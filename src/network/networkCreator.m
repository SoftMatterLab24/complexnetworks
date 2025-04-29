clear all
close all

%% 0. Settings

%%% Network properties
crystallinity = 1; % a value between 0 and 1

%%% Domain
L0 = 30;
xi = 1;
Rcut = 2*xi;

ClampWidth = 0.05; %percent of total height

DomainType = 'fixed'; %{'periodic','fixed'}

maxbond = 8;

filename = 'Network.dat';
%% I. Call the packing algorithm
% atoms: [ID | type | mol | x | y | z ]

[atoms,DomainBoundaries] = getAtoms(L0,xi,crystallinity,ClampWidth,DomainType);

%% II. Call the bonding algorithm

% 1. Construct the pairlist
[Idx] = buildPairlist(atoms,DomainBoundaries,Rcut);

% 2. Add bonds
[bonds] = getBonds(atoms,Idx,maxbond);

%% III. Write to .dat file
writeDat(filename,atoms,bonds,DomainBoundaries)
