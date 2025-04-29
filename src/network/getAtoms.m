function [atoms,DomainBoundaries] = getAtoms(L0,xi,crystallinity,ClampWidth,DomainType)

%%%%%%%%%%%%%%
%Density Based Approach
iplot = 1;

%%%%%%%%%%%%%%%

% User Parameters
rho = xi^2;                      % Number density

% Computed parameters
Ntarget = crystallinity*ceil(L0^2/rho);   % Target number of cross-link points
Ntarget
%Properties (these dont matter (remove in latter versions)
diap = 1;
rhop = 0.15;

%%% Size of domain (m) %%%
Lx  = 1*L0*xi;         % Length of x dimension
Ly  = 1*L0*xi;         % Length of y dimension

%%% Limits for x and y cps %%%
switch DomainType
    case 'periodic'
        ax = (Lx)/2;
        ay = (Ly)/2;
        bx = -ax;
        by = -ay;

        clamp_upperlimit = ay;
        clamp_lowerlimit = by;
    case 'fixed'
        ax = (Lx)/2;
        ay = (Ly)/2;
        bx = -ax;
        by = -ay;

        clamp_upperlimit = ay*(1 - ClampWidth);
        clamp_lowerlimit = by*(1 - ClampWidth);
end

%%% Random cp %%%
cp = [ax + (bx-ax)*rand(1) ay + (by-ay)*rand(1)];

%initialize
coords = cp;

itry = 1;
while 1 == 1

    if itry > 1e5
        break
    end

    cp = [ax + (bx-ax)*rand(1) ay + (by-ay)*rand(1)];

    % Get the distance btw randomly chosen cp and rest of cps
    data = repmat(cp',[1,size(coords,1)])'-coords;
    dist2 = data(:,1).^2+data(:,2).^2;

    % Get the indicies of any cps within one particle diameter
    iclose = find(dist2 < xi^2);

    % If there are no potential intersections - successful placement
    if  isempty(iclose)
        coords = [coords; cp];
        itry = 1;
        continue
    end

    itry = itry + 1;
    if size(coords,1) > Ntarget
        break
    end
end

natoms  = size(coords,1)

for ii = 1:natoms
    if any([coords(ii,2) > clamp_upperlimit, coords(ii,2) < clamp_lowerlimit])
        mol = 2;
    else
        mol = 1;
    end
      atoms(ii,:) = [ii 1 mol diap, rhop,coords(ii,1),coords(ii,2),0.0];
%     fprintf(fid,'%d %d %d %g %g %g %g %g\n',ii,ii,1,diap,rhop,coords(ii,1),coords(ii,2),0.0);
 end

if iplot
    figure(2); clf
    hold on
    scatter(coords(:,1),coords(:,2),36,'filled')
    axis equal
end

DomainBoundaries = [ax bx ay by];
% %%% Write to .dat file %%%
% disp('Writing to .dat file')
% 
% % Begin writing
% fid = fopen('Network1.dat','w');
% fprintf(fid,'Sacrifical Network\n\n');
% 
% % Atoms and bonds
% fprintf(fid,'%d atoms\n',natoms);
% % fprintf(fid,'%d bonds\n',nbonds);
% fprintf(fid,'2 atom types\n');
% fprintf(fid,'3 bond types\n');
% 
% % Simulation boundaries
% fprintf(fid,'%g %g xlo xhi\n',-Lx/2,Lx/2);
% fprintf(fid,'%g %g ylo yhi\n',-Ly/2,Ly/2);
% fprintf(fid,'%g %g zlo zhi\n\n',-Lx/20,Lx/20);
% 
% % For: atom_style hybrid bpm/sphere
% % atom-ID molecule-ID atom-type diameter density x y z
% fprintf(fid,'Atoms\n\n');
% for ii = 1:natoms
%     fprintf(fid,'%d %d %d %g %g %g %g %g\n',ii,ii,1,diap,rhop,coords(ii,1),coords(ii,2),0.0);
% end
% 
% natoms  = size(coords,1);
% rhotrue = natoms/Lx/Ly;
% 
% fprintf('Goal: %g, Achieved: %g\n',rho,rhotrue)

end
