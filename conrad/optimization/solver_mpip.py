"""
Define solver using the :mod:`cvxpy` module, if available.

For np.information on :mod:`cvxpy`, see:
http://www.cvxpy.org/en/latest/

If :func:`conrad.defs.module_installed` routine does not find the module
:mod:`cvxpy`, the variable ``SolverCVXPY`` is still defined in this
module's namespace as a lambda returning ``None`` with the same method
signature as the initializer for :class:`SolverCVXPY`. If :mod:`cvxpy`
is found, the class is defined normally.

Attributes:
	SOLVER_DEFAULT (:obj:`str`): Default solver, set to 'SCS' if module
		:mod:`scs` is installed, otherwise set to 'ECOS'.
"""
"""
Copyright 2016 Baris Ungun, Anqi Fu

This file is part of CONRAD.

CONRAD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CONRAD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CONRAD.  If not, see <http://www.gnu.org/licenses/>.
"""
from conrad.compat import *

import time
import numpy as np

from conrad.defs import vec as conrad_vec, module_installed, println
from conrad.medicine.dose import Constraint, MeanConstraint, MinConstraint, \
								 MaxConstraint, PercentileConstraint
from conrad.medicine.anatomy import Anatomy
from conrad.optimization.preprocessing import ObjectiveMethods
from conrad.optimization.solver_base import *

from numpy.random import rand


class SolverMPIP(Solver):
		"""
		Interface between :mod:`conrad` and :mod:`cvxpy` optimization library.

		:class:`SolverCVXPY` interprets :mod:`conrad` treatment planning
		problems (based on structures with attached objectives, dose
		constraints, and dose matrices) to build equivalent convex
		optimization problems using :mod:`cvxpy`'s syntax.

		The class provides an interface to modify, run, and retrieve
		solutions from optimization problems that can be executed on
		a CPU (or GPU, if :mod:`scs` installed with appropriate backend
		libraries).

		Attributes:
			problem (:class:`cvxpy.Minimize`): CVXPY representation of
				optimization problem.
			constraint_dual_vars (:obj:`dict`): Dictionary, keyed by
				constraint ID, of dual variables associated with each
				dose constraint in the CVXPY problem representation.
				The dual variables' values are stored here after each
				optimization run for access by clients of the
				:class:`SolverCVXPY` object.
		"""

		def __init__(self, n_beams=None, **options):
			"""
			Initialize empty :class:`SolverCVXPY` as :class:`Solver`.

			If number of beams provided, initialize the problem's
			:mod:`cvxpy` representation.

			Arguments:
				n_beams (:obj:`int`, optional): Number of beams in plan.
				**options: Arbitrary keyword arguments, passed to
					:meth:`SolverCVXPY.init_problem`.
			"""
			Solver.__init__(self)
			self.problem = None
			self.structures = None
			self.__x = None
			self.__constraint_indices = {}
			self.constraint_dual_vars = {}
			self.__solvetime = np.nan
			self.__status = None
			self.__value = None

			if isinstance(n_beams, int):
				self.init_problem(n_beams, **options)

		def init_problem(self, n_beams, use_slack=True, use_2pass=False,
						 **options):
			"""
			Initialize :mod:`cvxpy` variables and problem components.

			Create a :class:`cvxpy.Variable` of length-``n_beams`` to
			representthe beam  intensities. Invoke
			:meth:`SolverCVXPY.clear` to build minimal problem.

			Arguments:
				n_beams (:obj:`int`): Number of candidate beams in plan.
				use_slack (:obj:`bool`, optional): If ``True``, next
					invocation of :meth:`SolverCVXPY.build` will build
					dose constraints with slack variables.
				use_2pass (:obj:`bool`, optional): If ``True``, next
					invocation of :meth:`SolverCVXPY.build` will build
					percentile-type dose constraints as exact
					constraints instead of convex restrictions thereof,
					assuming other requirements are met.
				**options: Arbitrary keyword arguments.

			Returns:
				None
			"""
			self.__x = np.empty(n_beams, dtype=int)
			self.clear()

			self.use_slack = use_slack
			self.use_2pass = use_2pass
			self.gamma = options.pop('gamma', GAMMA_DEFAULT)

		@property
		def n_beams(self):
			""" Number of candidate beams in treatment plan. """
			return self.__x.size

		def clear(self):
			r"""
			Reset :mod:`cvxpy` problem to minimal representation.

			The minmal representation consists of:
				- An empty objective (Minimize 0),
				- A nonnegativity constraint on the vector of beam intensities (:math:`x \ge 0`).

			Reset dictionaries of:
				- Slack variables (all dose constraints),
				- Dual variables (all dose constraints), and
				- Slope variables for convex restrictions (percentile dose constraints).
			"""
			self.problem = None
			self.structures = None
			self.dvh_vars = {}
			self.slack_vars = {}
			self.constraint_dual_vars = {}
			self.__status = None
			self.__value = None

		@property
		def x(self):
			""" Vector variable of beam intensities, x. """
			return conrad_vec(self.__x)

		@property
		def x_dual(self):
			""" Dual variable corresponding to constraint x >= 0. """
			try:
				return conrad_vec(self.problem.constraints[0].dual_value)
			except:
				return None
		
		@property
		def x_var(self):
			return self.__x

		@property
		def solvetime(self):
			""" Solver run time. """
			return self.__solvetime

		@property
		def status(self):
			""" Solver status. """
			return self.__status

		@property
		def objective_value(self):
			""" Objective value at end of solve. """
			return self.__value

		@property
		def solveiters(self):
			""" Number of solver iterations performed. """
			return self.__solveiters

		def __objective_expression(self, structure):
			structure.normalize_objective()
			if structure.collapsable:
				return structure.objective.expr(structure.A_mean.T * self.__x)
			else:
				return structure.objective.expr(
						structure.A * self.x, structure.voxel_weights)

		def build(self, structures, exact=False, **options):
			"""
			Update :mod:`cvxpy` optimization based on structure data.

			Extract dose matrix, target doses, and objective weights
			from structures.

			Use doses and weights to add minimization terms to
			:attr:`SolverCVXPY.problem.objective`. Use dose constraints
			to extend :attr:`SolverCVXPY.problem.constraints`.

			(When constraints include slack variables, a penalty on each
			slack variable is added to the objective.)

			Arguments:
				structures: Iterable collection of :class:`Structure`
					objects.

			Returns:
				:obj:`str`: String documenting how data in
				``structures`` were parsed to form an optimization
				problem.
			"""
			self.clear()
			if isinstance(structures, Anatomy):
				structures = structures.list

			self.structures = structures


			# A, dose, weight_abs, weight_lin = \
					# self._Solver__gather_matrix_and_coefficients(structures)
			

			# self.problem.objective = cvxpy.Minimize(
			# 		weight_abs.T * cvxpy.abs(A * self.__x - dose) +
			# 		weight_lin.T * (A * self.__x - dose))

			# for s in structures:
				# self.__add_constraints(s, exact=exact)

			return self._Solver__construction_report(structures)

		def Qtilde(self, x, lbd):
			f = 0
			for s in self.structures:
				y = s.A @ x
				f = f + s.objective.eval(y)
			f = f + lbd * np.count_nonzero(x)
			return f

		def minSeed(self, set, x, ej, lbd):
			set = np.array(list(set))
			Qtilde = np.empty([set.size, 2])
			for j in range(set.size):
				e = np.zeros(self.n_beams)
				e[set[j]]= ej
				Qtilde[j,:] = [set[j], self.Qtilde(x + e, lbd)]
			jstar = np.argmin(Qtilde[:,1])
			return int(Qtilde[jstar,0]), Qtilde[jstar,1]

		def solve(self, **options):
			"""
			Execute optimization of a previously built planning problem.

			Arguments:
				**options: Keyword arguments specifying solver options,
					passed to :meth:`cvxpy.Problem.solve`.

			Returns:
				:obj:`bool`: ``True`` if :mod:`cvxpy` solver converged.

			Raises:
				ValueError: If specified solver is neither 'SCS' nor
					'ECOS'.
			"""

			# set verbosity level
			VERBOSE = bool(options.pop('verbose', VERBOSE_DEFAULT))
			PRINT = println if VERBOSE else lambda msg : None

			# solver options
			##solver = options.pop('solver', SOLVER_DEFAULT)
			##reltol = float(options.pop('reltol', RELTOL_DEFAULT))
			maxiter = int(options.pop('maxiter', MAXITER_DEFAULT))
			##use_gpu = bool(options.pop('gpu', GPU_DEFAULT))
			##use_indirect = bool(options.pop('use_indirect', INDIRECT_DEFAULT))

			# solve
			PRINT('running solver...')
			start = time.clock()
			x = np.zeros(self.n_beams)
			Shat = set(range(self.n_beams))
			S = set()
			lbd = 0.001
			Q0 = 0.02
			Qtilde = 1000.0
			Qtildeold = 1001.0
			k = 1
			while (Q0 < Qtilde) & (Qtilde < Qtildeold) & (k<maxiter):
			#for z in range(100):
				Qtildeold = Qtilde
				jstar, Qtilde = self.minSeed(Shat - S, x, 1, lbd)
				x[jstar] = x[jstar]	+ 1
				Sold = S.copy()
				S.add(jstar)
				if len(Sold) > 0:
					rstar, Qtilde2 = self.minSeed(Sold, x, -1, lbd)
					if Qtilde2 < Qtilde:
						x[rstar] = x[rstar] - 1
						S.remove(rstar)
						Qtilde = Qtilde2
				k = k+1
			self.__x = x
			self.__solvetime = time.clock() - start
			self.__solveiters = k
			self.__status = 'optimal'
			self.__value = self.Qtilde(x, lbd)


			PRINT("status: {}".format(self.__status))
			PRINT("optimal value: {}".format(self.__value))

			return True