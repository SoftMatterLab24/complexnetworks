function [bonds] = getBonds(atoms,Idx,maxbond)

natoms = length(atoms);


% for now assume random probability
P = 0.5;

fprintf('| Timestep:    | Bonds:     | \n')

bondsperatom = zeros([natoms,1]);
bonds = [];

N = 0; Ntry = 1000; nbonds = 0;
while 1 == 1
    N = N + 1;
    
    if N > Ntry
        fprintf('| Timestep: %4.0i | Bonds: %4.0i| \n',N,nbonds)
        break
    end

    %loop through atoms list
     for iatom = 1:length(atoms)
        
        % check to see if maxbond count is reached for atom ii
        if bondsperatom(iatom) >= maxbond-1
            continue
        end
        
        % see if any valid neighbors
        if isempty(Idx{iatom})
            continue
        end
        
        % draw random neighbor index
        ipick = randi(length(Idx{iatom}));
        jatom = Idx{iatom}(ipick);

        % ignore self
        if iatom == jatom
            continue
        end
        
        % check to see if maxbond count is reached for atom jj
        if bondsperatom(jatom) >= maxbond-1
            continue
        end

        % draw random probability 
        if rand(1) < P
            continue
        end
        
        %Add new bond to list
        bentry = [nbonds+1 1 iatom jatom]; %create bond entry
        bonds = [bonds; bentry]; %append 

        % Update bondcounters
        bondsperatom(iatom) = bondsperatom(iatom) + 1;
        bondsperatom(jatom) = bondsperatom(jatom) + 1;

        [nbonds, ~] = size(bonds);

        % Remove options from pair list
        Idx{iatom}(ipick) = [];
        Idx{jatom}(Idx{jatom} == iatom) = [];
        
     end

        fprintf('| Timestep: %4.0i | Bonds: %4.0i | Max BpA: %4.0i | Mean BpA %4.4f |\n',N,nbonds,max(bondsperatom),mean(bondsperatom))
        %pause(0.1)
 end
end