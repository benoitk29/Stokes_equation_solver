function [tilde_AA, tilde_LL] = elimine_stokes(AA, LL, Refneu, Coorneu,probleme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Réalisation de la pseudo-élimination des conditions aux limites dans le
% système d'équations de Stokes
%
% SYNOPSIS [tilde_AA, tilde_LL] = elimine_stokes(AA, LL, Refneu, Coorneu)
%
% INPUT
% -----
% * AA : matrice du système linéaire AU=L (dimension (2N+Ns)x(2N+Ns))
%
% * LL : vecteur des termes sources (dimension (2N+Ns)x1)
%
% * Refneu : vecteur permettant de reconnaitre les noeuds du bord et quelle
% condition appliquée
%
% * Coorneu : matrice contenant les coordonnées (x, y) de chaque point
%   du maillage.
%
% * probleme  : type de problème considéré ('rectangle' ou 'marche').
%
% OUTPUT
% ------
% - tilde_AA : matrice du système linéaire modifiée, après application des conditions
%   aux limites.
%
% - tilde_LL : vecteur des termes sources modifiée, après application des
%   conditions aux limites.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nombre_points, ~] = size(Refneu);

% On initialise les matrices tilde_AA et tilde_LL qui contiendront les matrices modifiées
tilde_AA = AA;
tilde_LL = LL;

% Boucle sur chaque point pour traiter les conditions aux limites
for idx = 1:nombre_points
    x=Coorneu(idx, 1);
    y=Coorneu(idx, 2);
    if strcmp(probleme, 'rectangle')
        % Si refneu égal à 1, cela signifie qu'on est sur une condition
        % de Dirichlet, où les vitesses doivent être fixées (u = g).
        if Refneu(idx) == 1
            % Annulation des lignes et colonnes de la matrice AA correspondantes au point
            % pour que la vitesse soit imposée.
            tilde_AA(idx, :) = 0;
            tilde_AA(idx + nombre_points, :) = 0;

            % Pour  que la vitesse reste fixée aux valeurs (u1 et u2),on met 1
            % sur les termes diagonaux correspondants (cf amphi 5) (afin
            % d'avoir u1=g1, ...)
            tilde_AA(idx, idx) = 1;
            tilde_AA(idx + nombre_points, idx + nombre_points) = 1;

            % Application des valeurs des conditions aux limites de Dirichlet pour les
            % vitesses u1 et u2 (en fonction des coordonnées du point).
            tilde_LL(idx) = g1(x,y,probleme); % Pour u1
            tilde_LL(idx + nombre_points) = g2(x,y); % Pour u2
        end

        %  Si refneu égal à 2, cela signifie qu'on est sur une condition
        % de Dirichlet
        if Refneu(idx) == 2
            % Annulation des lignes et colonnes de la matrice AA associées à ce point,
            % car ces points n'ont pas d'inconnue associée, seulement une condition de Dirichlet.
            tilde_AA(idx, :) = 0;
            tilde_AA(:, idx) = 0;
            tilde_AA(idx + nombre_points, :) = 0;
            tilde_AA(:, idx + nombre_points) = 0;

            % Pour  que la vitesse reste fixée aux valeurs (u1 et u2),on met 1
            % sur les termes diagonaux correspondants (cf amphi 5) (afin
            % d'imposer la condition à la limite)
            tilde_AA(idx, idx) = 1;
            tilde_AA(idx + nombre_points, idx + nombre_points) = 1;

            % On met à zéro les valeurs du vecteur LL associées à ce point,
            % pour des conditions de Dirichlet homogènes
            tilde_LL(idx) = 0;            % Pour u1
            tilde_LL(idx + nombre_points) = 0;  % Pour u2
        end



    elseif strcmp(probleme, 'marche')

        % Si Refneu égal à 1 ou 3, cela signifie qu'on impose une condition
        % de Dirichlet où les vitesses doivent être nulles (Γ_up ou Γ_down).
        if ((Refneu(idx) == 1) || (Refneu(idx) == 3))
            % Annulation des lignes et colonnes correspondantes dans la matrice
            tilde_AA(idx, :) = 0;
            tilde_AA(idx + nombre_points, :) = 0;

            % On impose 1 sur les termes diagonaux pour fixer les vitesses à 0
            tilde_AA(idx, idx) = 1;
            tilde_AA(idx + nombre_points, idx + nombre_points) = 1;

            % Les conditions de Dirichlet imposent u1 = 0 et u2 = 0
            tilde_LL(idx) = 0; % Pour u1
            tilde_LL(idx + nombre_points) = 0; % Pour u2

            % Si Refneu égal à 2, on est sur Γ_right avec une condition de Neumann homogène.
            % La condition sur le bord droit est naturellement imposée par
            % la formule variationnelle
        elseif Refneu(idx) == 4
            % Annulation des lignes et colonnes correspondantes dans la matrice
            tilde_AA(idx, :) = 0;
            tilde_AA(idx + nombre_points, :) = 0;

            % On impose 1 sur les termes diagonaux pour fixer les valeurs
            tilde_AA(idx, idx) = 1;
            tilde_AA(idx + nombre_points, idx + nombre_points) = 1;

            % Conditions imposées via les fonctions g1 et g2
            tilde_LL(idx) = g1(x,y, probleme); % Pour u1
            tilde_LL(idx + nombre_points) = g2(x,y); % Pour u2

        end
    else
        error('Type de problème non reconnu par la fonction g1. Utilisez ''rectangle'' ou ''marche''.');
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%24

