# MIP_MCLPD
## Notation:
RS: Recharge stations \
DC: Distribution centers\
Facility: Union of recharge stations and distribution centers 
## Indices:
i,k: location indices for distribution centers and recharge stations\
j: demand index, j=1,...,m
## Parameters:
f<sub>d</sub>: Maximum flight range for delivering demand and return to the facility\
f<sub>p</sub>: Maximum flight range for reaching other facilities safely\
a<sub>j</sub> : Rate of demand at demand j\
c<sub>DC</sub>: cost of opening a distribution center\
c<sub>RS</sub>: cost of opening a recharge station
## Sets:
P<sub>i</sub>: Set of facility locations reachable from facility location i (distance between facilities <= fp) \
R<sub>j</sub> : Set of facility locations reachable from demand j (distance between facility and demand <= fd) \
L<sub>RS</sub>: Set of possible locations for recharge stations \
L<sub>DC</sub>: Set of possible locations for distribution centers

![math_formula](https://user-images.githubusercontent.com/53580699/134614272-5496fdd6-2705-4a3b-97f9-1755d19eff51.png)

Objective *(1)* maximizes total satisfed demand. Objective *(2)* minimizes
total cost of opening distribution centers and recharging stations. Con-
straint *(3)* ensures that demands can only be satisfed from opened facil-
ity which is reachable. Constraint *(4)* provides there should be outflow if a recharge station is placed. Constraint *(5)* guarantees that there is inflow
if a distribution center is placed. Constraint *(6)* ensures that outflow from
a recharge station should be greater than inflow. Set *(7)* represent binary
location decision and *(8)* guarantees that fows are integer.
