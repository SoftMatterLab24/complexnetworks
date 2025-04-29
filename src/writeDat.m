function [] = writeDat(filename,atoms,bonds,DomainBoundaries)

[natoms ~] = size(atoms);
[nbonds ~] = size(bonds);

diap = 0.5;
rhop = 1;

%%% Write to .dat file %%%
disp('Writing to .dat file')

% Begin writing
fid = fopen(filename,'w');
fprintf(fid,'Sacrifical Network\n\n');

% Atoms and bonds
fprintf(fid,'%d atoms\n',natoms);
fprintf(fid,'%d bonds\n',nbonds);
fprintf(fid,'%d atom types\n',max(atoms(:,2)));
fprintf(fid,'1 bond types\n');

% Simulation boundaries
fprintf(fid,'%g %g xlo xhi\n',DomainBoundaries(2),DomainBoundaries(1));
fprintf(fid,'%g %g ylo yhi\n',DomainBoundaries(4),DomainBoundaries(3));
fprintf(fid,'%g %g zlo zhi\n\n',DomainBoundaries(2)/20,DomainBoundaries(1)/20);

% For: atom_style hybrid bpm/sphere
% atom-ID molecule-ID atom-type diameter density x y z
fprintf(fid,'Atoms\n\n');
for ii = 1:natoms
    fprintf(fid,'%d %d %d %g %g %g %g %g\n',atoms(ii,1),atoms(ii,3),atoms(ii,2),diap,rhop,atoms(ii,6),atoms(ii,7),atoms(ii,8));
end

fprintf(fid,'\nBonds\n\n');
for jj = 1:nbonds
    fprintf(fid,'%d %d %d %d\n',bonds(jj,1),bonds(jj,2),bonds(jj,3),bonds(jj,4));
end

fclose(fid);
disp('done')
%%
%%% VISUALIZTION %%% 
%%% Write to .dat file %%%
disp('Writing to visualization file')

% Begin writing
fid = fopen('NetworkVisual.dat','w');
fprintf(fid,'Sacrifical Network\n\n');

% Atoms and bonds
fprintf(fid,'%d atoms\n',natoms);
fprintf(fid,'%d bonds\n',nbonds);
fprintf(fid,'%d atom types\n',max(atoms(:,2)));
fprintf(fid,'1 bond types\n');

% Simulation boundaries
fprintf(fid,'%g %g xlo xhi\n',DomainBoundaries(2),DomainBoundaries(1));
fprintf(fid,'%g %g ylo yhi\n',DomainBoundaries(4),DomainBoundaries(3));
fprintf(fid,'%g %g zlo zhi\n\n',DomainBoundaries(2)/20,DomainBoundaries(1)/20);

% For: atom_style hybrid bpm/sphere
% atom-ID molecule-ID atom-type diameter density x y z
fprintf(fid,'Atoms\n\n');
for ii = 1:natoms
    fprintf(fid,'%d %d %g %g %g %d %g %g \n',atoms(ii,1),atoms(ii,2),atoms(ii,6),atoms(ii,7),atoms(ii,8),atoms(ii,3),diap,rhop);
end

fprintf(fid,'\nBonds\n\n');
for jj = 1:nbonds
    fprintf(fid,'%d %d %d %d\n',bonds(jj,1),bonds(jj,2),bonds(jj,3),bonds(jj,4));
end

fclose(fid);
disp('done')
end