# =========================
# DMGB.mod   (AMPL model)
# =========================

set G;                 # generators
set T0 ordered;         # time periods including 0
set T within T0;        # time periods excluding 0
set S;                 # scenarios

param pi {S} > 0;       # scenario probabilities
param lambda {T,S};     # day-ahead prices (EUR/MWh)

param cq {G} >= 0;
param cl {G} >= 0;
param cb {G} >= 0;

param Pmin {G} >= 0;
param Pmax {G} >= 0;

param cU {G} >= 0;
param cD {G} >= 0;

param tU {G} integer >= 1;   # minimum up time
param tD {G} integer >= 1;   # minimum down time

param u0 {G} binary;         # initial commitment at t=0

var u {G,T0} binary;
var v {G,T0} binary;         # start-up
var w {G,T0} binary;         # shut-down

var P {G,T,S} >= 0;          # dispatched power (MW)

# Minimize expected negative profit (convex MIQP form)
minimize Expected_Neg_Profit:
    sum {t in T, s in S} pi[s] *
        sum {i in G} (
            cq[i]*P[i,t,s]^2 + cl[i]*P[i,t,s] + cb[i]*u[i,t]
            - lambda[t,s]*P[i,t,s]
        )
  + sum {t in T0} sum {i in G} ( cU[i]*v[i,t] + cD[i]*w[i,t] );

subject to Initial_Commitment {i in G}:
    u[i,0] = u0[i];

subject to Commitment_Logic {i in G, t in T}:
    u[i,t] - u[i,prev(t)] - v[i,t] + w[i,t] = 0;

# Minimum up-time (standard window constraint)
subject to Min_Up_Time {i in G, t in T}:
    sum {k in T: ord(k) <= ord(t) and ord(k) >= ord(t) - tU[i] + 1} v[i,k] <= u[i,t];

# Minimum down-time (standard window constraint)
subject to Min_Down_Time {i in G, t in T}:
    sum {k in T: ord(k) <= ord(t) and ord(k) >= ord(t) - tD[i] + 1} w[i,k] <= 1 - u[i,t];

subject to Min_Generation {i in G, t in T, s in S}:
    P[i,t,s] >= Pmin[i]*u[i,t];

subject to Max_Generation {i in G, t in T, s in S}:
    P[i,t,s] <= Pmax[i]*u[i,t];
