# CRTBP-Library

This code serves as a matlab Library for non-dimensional CRTBP System. Obtains: Lyapunov Orbit and Halo Orbit Parameters. Future Updates: Functions for all symmetric and asymmetric orbits (given in Daniel Grebow's Master's thesis) , Invarient Manifold functions, Resonant Orbit Functions (as in Manninie Gupta's Master's Thesis) and more. Invarient Manifold Code can be found in 3-Body_Problem Code.


## Functions

Make use of following functions:

* `LyapOrbitParameters.m` gives the Initial Conditions and other parameters (time period, Eigen Structure, Stability Index etc.) in the non-dimentinal CR3BP of a lyapunov orbit for spacific mu, liberation point and jacobian.
* `LyapOrbitFamilyParameters.m` gives the Initial Conditions and other parameters (time period, Eigen Structure, Stability Index etc.) in the non-dimentinal CR3BP of family of lyapunov orbits for spacific mu and liberation point.
* `HaloOrbitParameters.m` gives the Initial Conditions and other parameters (time period, Eigen Structure, Stability Index etc.) in the non-dimentinal CR3BP of a Halo orbit (northern and Southern should be specified) for spacific mu, liberation point and jacobian.
* `HaloOrbitFamilyParameters.m` gives the Initial Conditions and other parameters (time period, Eigen Structure, Stability Index etc.) in the non-dimentinal CR3BP of family of Halo orbits (northern and Southern should be specified) for spacific mu and liberation point.
* `PeriodicOrbitInvariantMfdsIC.m` returns the Initial Conditions to compute Invariant Manifolds (Stable+/-, Unstable+/-) for a given Periodic Orbit. The direction needs to be taken care of!
* `Integrator.m` make use of this fucntion to compute the final Periodic Orbits

## Examples

Inputting System Variables

```ruby 
UserDat.Dimension      = 3;
UserDat.mu        = 0.0121505856;
UserDat.PointLoc       = 1;
UserDat.NoOfFam        = 50;
UserDat.CorrectionPlot = 1;
G_var                  = GlobalData(UserDat);
```

- `UserDat.PointLoc` specifies which equilibrium points you want the data for? can take 1,2,3 for now
- `UserDat.CorrectionPlot` can take 1 and 0 specifying if you want differential correction plots or not.

```ruby
e = 3.16;
[LyapOrb] = LyapOrbitParameters(UserDat,G_var,e);
figure()
[t,x] = Integrator(G_var,G_var.IntFunc.EOM,LyapOrb.IC,[0 LyapOrb.time],'forward');
plot(x(:,1),x(:,2));hold on; grid on;
scatter(1-UserDat.mu,0,'p','filled');
xlabel('x (ND)')
ylabel('y (ND)')
title('Lyapunov Orbit Family at L1 for mu = 0.01215 and e = 3.16');
```

```ruby
e = 3.16;
[LyapOrb] = LyapOrbitParameters(UserDat,G_var,e);
figure()
[t,x] = Integrator(G_var,G_var.IntFunc.EOM,LyapOrb.IC,[0 LyapOrb.time],'forward');
plot(x(:,1),x(:,2));hold on; grid on;
scatter(1-UserDat.mu,0,'p','filled');
xlabel('x (ND)')
ylabel('y (ND)')
title('Lyapunov Orbit Family at L1 for mu = 0.01215 and e = 3.16');
```


### Example:1 Plotting Lyapunov Orbit Family


```ruby
[LyapOrbFam] = LyapOrbitFamilyParameters(UserDat,G_var);
figure()
plot(LyapOrbFam.Energy, LyapOrbFam.time,'-r','LineWidth',2);
figure()
set(0,'DefaultAxesColorOrder',flipud(jet(length(LyapOrbFam.Energy))));
for i = 1:length(LyapOrbFam.Energy)
    [t,x] = Integrator(G_var,G_var.IntFunc.EOM,LyapOrbFam.IC(i,:),[0 LyapOrbFam.time(i)],'forward');
    plot(x(:,1),x(:,2));hold on; grid on;
end
scatter(1-UserDat.mu,0,'p','filled');
caxis([LyapOrbFam.Energy(end), LyapOrbFam.Energy(1)]);colormap("jet");
colorbar
xlabel('x (ND)')
ylabel('y (ND)')
title('Lyapunov Orbit Family at L1 for mu = 0.01215');

```
<img src="images/LyapOrbit1.png" width="400">
You can also carry out the Jacobian and Time Period Analysis
<img src="images/LyapOrbitPeriodStudy.png" width="400">
Similarly you can carry out a Stability Analysis using Stability Index.

### Example:2 Plotting Halo Orbit Family

```ruby
[HaloOrbFam]              = HaloOrbitFamilyParameters(UserDat,G_var, 'northern');
figure()
set(0,'DefaultAxesColorOrder',flipud(jet(length(HaloOrbFam.Energy))));
for i = 1:length(HaloOrbFam.Energy)
    [t,x] = Integrator(G_var,G_var.IntFunc.EOM,HaloOrbFam.IC(i,:),[0 HaloOrbFam.time(i)],'forward');
    plot3(x(:,1),x(:,2),x(:,3));hold on; grid on;
end
scatter3(1-UserDat.mu,0,0,'p','filled');
caxis([HaloOrbFam.Energy(end), HaloOrbFam.Energy(1)]);
colorbar
xlabel('x (ND)')
ylabel('y (ND)')
zlabel('z (ND)')
title('Lyapunov Orbit Family at L1 for mu = 0.01215');
```
<img src="images/northernHalo1.png" width="400">
<img src="images/northernHalo2.png" width="400">

### Example:3 Plotting Invariant Manifolds

```ruby
[LyapOrb]              = LyapOrbitParameters(UserDat,G_var, 3.17);
ans1 = PeriodicOrbitInvariantMfdsIC(G_var,LyapOrb,80, 'stable',1);
figure()
for i = 1:length(ans1(:,1))
    [t,x] = Integrator(G_var,G_var.IntFunc.EOM,ans1(i,1:6),[0 1.5*LyapOrb.time],'backward');
    plot(x(:,1),x(:,2),'g');hold on; grid on;
end
scatter(1-UserDat.mu,0,'p','filled','r');
xlabel('x (ND)')
 ylabel('y (ND)')
title('Invariant Manifold for L1 Lyapunov Orbit');
subtitle('mu = 0.01215, Jacobian=3.17');
```
<img src="images/LyapOrbitInvariantManifold.png" width="400">


## References Used


    - Dynamical Syatems, the three-body problem and Space Mission Design, Koon, Lo, Marsden, Ross
    - Generating Periodic Orbits In CR3BP With Applications To Lunar South Pole Coverage - Daniel Grebow
    - Finding Order in Chaos: Resonant Orbits and Poincare Section - Maaninee Gupta
    - TADPOLE ORBITS IN THE L4/L5 REGION: CONSTRUCTION AND LINKS TO OTHER FAMILIES OF PERIODIC ORBITS - Alexandre G. Van Anderlecht
    - SPACECRAFT TRAJECTORY DESIGN TECHNIQUES USING RESONANT ORBITS - Srianish Vutukuri
    - REPRESENTATIONS OF INVARIANT MANIFOLDS FOR APPLICATIONS IN SYSTEM-TO-SYSTEM TRANSFER DESIGN - Christopher E. Patterson

## More to Come...

More to be added :

* Halo Orbit for L3
* Vertical Orbit (L1/2/3/4/5)
* Axial Orbit (L1/2/3/4/5)
* Butterfly (L2/?)
* Planar (L4/5)
* Tadpole Orbit, Horseshoe Orbit (L4/5)
* Resonant Orbits
* Bifurcations and More...
