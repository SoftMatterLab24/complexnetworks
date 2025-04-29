function  [Idx] = buildPairlist(atoms,DomainBoundaries,Rcut)
    
    PairlistMethod = 1;

    CoordsA = [atoms(:,6) atoms(:,7)];
    CoordsB = [atoms(:,6) atoms(:,7)];

    [a,b] = size(CoordsB);


    % needed if periodic BC (flat)
    %if PairlistMethod == 2
    Ldom = DomainBoundaries(2)-DomainBoundaries(1);
    for ii = 1:length(CoordsA) %loop through lattice coords

            xi = CoordsA(ii,1); 
            yi = CoordsA(ii,2);
    
            intvector = zeros(length(CoordsA),1);
 
            for jj = 1:a %loop through robot coords
        
                    xj = CoordsB(jj,1); 
                    yj = CoordsB(jj,2);

                    dx = abs(xi-xj);
                    dy = abs(yi-yj);

                    if PairlistMethod == 2
                        if dx > Ldom/2
                            dx = Ldom-dx;
                        end
                        if dy > Ldom/2
                            dy = Ldom-dy;
                        end
                    end
                    D(ii,jj) = norm([dx dy]);

                    if D(ii,jj) < Rcut
                        intvector(jj) = true;
                    end

            end
    
            Idx{ii} = find(intvector==true);
    end

    Idx = Idx';

end