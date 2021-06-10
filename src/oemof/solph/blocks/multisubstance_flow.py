# -*- coding: utf-8 -*-

"""Creating sets, variables, constraints and parts of the objective function
for Flow objects.

SPDX-FileCopyrightText: Uwe Krien <krien@uni-bremen.de>
SPDX-FileCopyrightText: Simon Hilpert
SPDX-FileCopyrightText: Cord Kaldemeyer
SPDX-FileCopyrightText: Patrik Schönfeldt
SPDX-FileCopyrightText: Birgit Schachler
SPDX-FileCopyrightText: jnnr
SPDX-FileCopyrightText: jmloenneberga

SPDX-License-Identifier: MIT

"""

from pyomo.core import Constraint
from pyomo.core import NonNegativeReals
from pyomo.core import Set
from pyomo.core import Var
from pyomo.core.base.block import SimpleBlock


class MultiSubstanceFlow(SimpleBlock):
    r""" Flow block with definitions for standard flows.

    **The following variables are created**:

    negative_gradient :
        Difference of a flow in consecutive timesteps if flow is reduced
        indexed by NEGATIVE_GRADIENT_FLOWS, TIMESTEPS.

    positive_gradient :
        Difference of a flow in consecutive timesteps if flow is increased
        indexed by NEGATIVE_GRADIENT_FLOWS, TIMESTEPS.

    **The following sets are created:** (-> see basic sets at :class:`.Model` )

    SUMMED_MAX_FLOWS
        A set of flows with the attribute :attr:`summed_max` being not None.
    SUMMED_MIN_FLOWS
        A set of flows with the attribute :attr:`summed_min` being not None.
    NEGATIVE_GRADIENT_FLOWS
        A set of flows with the attribute :attr:`negative_gradient` being not
        None.
    POSITIVE_GRADIENT_FLOWS
        A set of flows with the attribute :attr:`positive_gradient` being not
        None
    INTEGER_FLOWS
        A set of flows where the attribute :attr:`integer` is True (forces flow
        to only take integer values)

    **The following constraints are build:**

    Flow max sum :attr:`om.Flow.summed_max[i, o]`
      .. math::
        \sum_t flow(i, o, t) \cdot \tau
            \leq summed\_max(i, o) \cdot nominal\_value(i, o), \\
        \forall (i, o) \in \textrm{SUMMED\_MAX\_FLOWS}.

    Flow min sum :attr:`om.Flow.summed_min[i, o]`
      .. math::
        \sum_t flow(i, o, t) \cdot \tau
            \geq summed\_min(i, o) \cdot nominal\_value(i, o), \\
        \forall (i, o) \in \textrm{SUMMED\_MIN\_FLOWS}.

    Negative gradient constraint
      :attr:`om.Flow.negative_gradient_constr[i, o]`:
        .. math::
          flow(i, o, t-1) - flow(i, o, t) \geq \
          negative\_gradient(i, o, t), \\
          \forall (i, o) \in \textrm{NEGATIVE\_GRADIENT\_FLOWS}, \\
          \forall t \in \textrm{TIMESTEPS}.

    Positive gradient constraint
      :attr:`om.Flow.positive_gradient_constr[i, o]`:
        .. math:: flow(i, o, t) - flow(i, o, t-1) \geq \
          positive\__gradient(i, o, t), \\
          \forall (i, o) \in \textrm{POSITIVE\_GRADIENT\_FLOWS}, \\
          \forall t \in \textrm{TIMESTEPS}.

    **The following parts of the objective function are created:**

    If :attr:`variable_costs` are set by the user:
      .. math::
          \sum_{(i,o)} \sum_t flow(i, o, t) \cdot variable\_costs(i, o, t)

    The expression can be accessed by :attr:`om.Flow.variable_costs` and
    their value after optimization by :meth:`om.Flow.variable_costs()` .

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _create(self, group=None):
        r"""Creates sets, variables and constraints for all standard flows.

        Parameters
        ----------
        group : list
            List containing tuples containing flow (f) objects and the
            associated source (s) and target (t)
            of flow e.g. groups=[(s1, t1, f1), (s2, t2, f2),..]
        """
        if group is None:
            return None

        m = self.parent_block()

        # ########################## SETS #####################################
        self.MULTISUBSTANCEFLOWS = Set(
            initialize=[(g[0], g[1]) for g in group])

        # ######################### Variables  ################################

        self.substance_flow = Var(
            m.MULTISUBSTANCEFLOWS, m.SUBSTANCES, m.TIMESTEPS,
            within=NonNegativeReals)

        # ######################### CONSTRAINTS ###############################

        def _substance_flow_rule(block, i, o, s, t):
            """Rule definition for multisubstances flows"""
            lhs = m.flow[i, o, t] * m.substances[s][t]
            rhs = block.substance_flow[i, o, s, t]
            return lhs == rhs

        self.substance_flow_rule = Constraint(
            m.MULTISUBSTANCEFLOWS, m.SUBSTANCES, m.TIMESTEPS,
            rule=_substance_flow_rule)
