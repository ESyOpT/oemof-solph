# -*- coding: utf-8 -*-

"""
solph version of oemof.network.energy_system

SPDX-FileCopyrightText: Uwe Krien <krien@uni-bremen.de>
SPDX-FileCopyrightText: Simon Hilpert
SPDX-FileCopyrightText: Cord Kaldemeyer
SPDX-FileCopyrightText: Stephan Günther
SPDX-FileCopyrightText: Birgit Schachler

SPDX-License-Identifier: MIT

"""

from oemof.network import energy_system as es


class EnergySystem(es.EnergySystem):
    """A variant of :class:`EnergySystem
    <oemof.core.energy_system.EnergySystem>` specially tailored to solph.

    In order to work in tandem with solph, instances of this class always use
    `solph.GROUPINGS <oemof.solph.GROUPINGS>`. If custom groupings are
    supplied via the `groupings` keyword argument, `solph.GROUPINGS
    <oemof.solph.GROUPINGS>` is prepended to those.

    If you know what you are doing and want to use solph without
    `solph.GROUPINGS <oemof.solph.GROUPINGS>`, you can just use
    :class:`core's EnergySystem <oemof.core.energy_system.EnergySystem>`
    directly.
    """

    def __init__(self, **kwargs):
        # Doing imports at runtime is generally frowned upon, but should work
        # for now. See the TODO in :func:`constraint_grouping
        # <oemof.solph.groupings.constraint_grouping>` for more information.
        from oemof.solph.groupings import GROUPINGS

        kwargs["groupings"] = GROUPINGS + kwargs.get("groupings", [])

        super().__init__(**kwargs)

        self.substances = kwargs.get('substances', set())

    def add(self, *nodes):
        """Add :class:`nodes <oemof.network.Node>` to this energy system."""
        # check substances for all added nodes to ensure only substances
        # declared in energySystem are used
        for n in nodes:
            # check flows for both inputs and outputs
            for f in (set(n.inputs.values()) | set(n.outputs.values())):
                if hasattr(f, 'substances') and f.substances:
                    if not set(f.substances.keys()).issubset(self.substances):
                        # if flow contains additional substances not
                        # defined in energysystem calculate set difference
                        # and inform user
                        s = set(f.substances.keys()) - self.substances
                        msg = ("Substance {} is  used in flow but not"
                               + " declared in solph.EnergySystem")
                        raise ValueError(msg.format(s))

        # add nodes if substances are declared correctly
        self.nodes.extend(nodes)
