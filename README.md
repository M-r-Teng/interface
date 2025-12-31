Calculating analytical and numerical solutions for FRP interfacal problems:

ana_from_Wu：The analytical solution proposed in Wu's paper, as shown in Figure 4 of this paper. [18] Wu YF, Xu XS, Sun JB, Jiang C. Analytical solution for the bond strength of externally bonded reinforcement. Composite Structures. 2012;94:3232-9. 

ana_linear：The analytical solution of the classical linear model shown in Figure 2 (a).

There is currently no analytical solution for Popovics model.

frpodes0con_popo: s0 control, Popovics bond-slip model
frpbvpslcon_popo: sl control, Popovics bond-slip model

frpodes0con_Wu: s0 control,  model in [18].
frpbvpslcon_Wu:sl control,  model in [18].

frpodes0con_linear：s0 control,  linear model.

HugoFcon：reproduced the force control algorithm based on Hugo C. Biscaia's multiple papers.

frpbvps0con_EC:  interface with lateral confinement, not shown in the paper. [58] Biscaia HC, Chastre C, Silva MAG. Bond-slip model for FRP-to-concrete bonded joints under external 
compression. Composites Part B: Engineering. 2015;80:246-59.

s_**： single slip distributions, including the force control condition.
