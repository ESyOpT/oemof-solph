# -*- coding: utf-8 -*-

"""Creating sets, variables, constraints and parts of the objective function
for Bus objects.

SPDX-FileCopyrightText: Uwe Krien <krien@uni-bremen.de>
SPDX-FileCopyrightText: Simon Hilpert
SPDX-FileCopyrightText: Cord Kaldemeyer
SPDX-FileCopyrightText: Patrik Schönfeldt
SPDX-FileCopyrightText: Birgit Schachler
SPDX-FileCopyrightText: jnnr
SPDX-FileCopyrightText: jmloenneberga

SPDX-License-Identifier: MIT

"""

from pyomo.core import BuildAction
from pyomo.core import Constraint
from pyomo.core.base.block import SimpleBlock
from oemof.tools import debugging
import warnings


class MultiSubstanceBus(SimpleBlock):
    r"""Block for all balanced buses.

    **The following constraints are build:**

    Bus balance  :attr:`om.Bus.balance[i, o, t]`
      .. math::
        \sum_{i \in INPUTS(n)} flow(i, n, t) =
        \sum_{o \in OUTPUTS(n)} flow(n, o, t), \\
        \forall n \in \textrm{BUSES},
        \forall t \in \textrm{TIMESTEPS}.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _create(self, group=None):
        """Creates the balance constraints for the class:`Bus` block.

        Parameters
        ----------
        group : list
            List of oemof bus (b) object for which the bus balance is created
            e.g. group = [b1, b2, b3, .....]
        """
        if group is None:
            return None

        m = self.parent_block()

        ins = {}
        outs = {}
        for n in group:
            ins[n] = [i for i in n.inputs]
            outs[n] = [o for o in n.outputs]
            # get the concentrations of first output flow as reference
            ref_conc = list(n.outputs.values())[0].substances.items()
            # are all output flows concentrations equal to ref_conc?
            result = all(elem == ref_conc for elem in [x.substances.items() for x in n.outputs.values()])
            if result:
                pass
            else:
                msg = "\n\nAll flows connected to the outputs of a " \
                        + "MultiSubstanceBus must have the same"\
                        + " substance concentrations. Results may otherwise "\
                        + "be wrong.\n"
                # get some infos for the warning message
                source = n.label
                output_labels = [x.label for x in n.outputs.keys()]
                concentrations = [dict(x.substances.items()) for x in n.outputs.values()]
                # append the warning message with those infos
                for i in range(len(output_labels)):
                    msg += "\nFlow: {} -> {}: {}".format(source,
                                                output_labels[i],
                                                concentrations[i])
                warnings.warn(msg, debugging.SuspiciousUsageWarning)

        def _busbalance_rule(block):
            for t in m.TIMESTEPS:
                for g in group:
                    for s in m.SUBSTANCES:
                        lhs = sum(
                            m.substance_flow[i, g, s, t] for i in ins[g])
                        rhs = sum(
                            m.substance_flow[g, o, s, t] for o in outs[g])
                        expr = lhs == rhs
                        # no inflows no outflows yield: 0 == 0 which is True
                        """
                        expr is a EqualityExpression object if lhs and rhs
                        are objects. The if clause does not check if the
                        value inside the Equality Expression is True, but if
                        the expr object is of type True. If it is of type
                        True, no inputs and outputs are present and no
                        constraint needs to be build. If it's not of type True,
                        a constraint needs to be build that ensures that the
                        sum of all input flows of a specific substance equals
                        the sum of all output flows of that substance.
                        """
                        if expr is not True:
                            # The following line adds a constraint with the
                            # current g, s and t with the rule lhs==rhs to the
                            # optimization
                            block.balance.add((g, s, t), expr)

        self.balance = Constraint(
            group, m.SUBSTANCES, m.TIMESTEPS, noruleinit=True)
        self.balance_build = BuildAction(rule=_busbalance_rule)
