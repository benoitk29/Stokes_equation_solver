function res = g1(x, y,probleme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la condition de Dirichlet pour la composante u1
%
% SYNOPSIS [res] = g1(x, y)
%
% INPUT * x, y : coordonnées du point dans le plan
%       * probleme  : type de problème considéré ('rectangle' ou 'marche').
%
% OUTPUT - res : valeur de la condition de Dirichlet pour u1,
%                en fonction des coordonées et du problème
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(probleme, 'rectangle')
    if x == 0
        res =(2 - y)*y;
    else
        res = 0;
    end
elseif strcmp(probleme, 'marche')

    res =  (y - 1) * (2 - y)/(0.25);
else
    error('Type de problème non reconnu par la fonction g1. Utilisez ''rectangle'' ou ''marche''.');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%24