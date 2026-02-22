function [Kel] = matK_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de rigidité elementaire en P2 Lagrange
%
% SYNOPSIS [Kel] = matK_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de rigidité elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de rigidité
% -------------------------------
% Initialisation
Kel =zeros(6, 6);  % A COMPLETER

% Points et poids de quadrature
S_hat = [1/6, 1/6;
    2/3, 1/6;
    1/6, 2/3];
poids = 1/6;

%Traitement de Bl pour calculer la transposée de son inverse utilisée dans
% le calcul de quadrature, tout comme son déterminant
B = [x2 - x1, x3 - x1;
    y2 - y1, y3 - y1];
invB = inv(B);
detB = abs(det(B));

for q = 1:size(S_hat, 1)
    xsi = S_hat(q, 1);
    eta = S_hat(q, 2);
    %Après avoir calculé chaque gradient de chaque fonction de bases wi, on
    %fait des 6 gradients trouvés une matrice 6x6 afin de calcul la matrice
    %élémentaire via un produit de matrices
    gradN_hat = [-3 + 4*xsi + 4*eta, 4*xsi - 1, 0, 4*(1 - 2*xsi - eta), 4*eta, -4*eta;
        -3 + 4*xsi + 4*eta, 0, 4*eta - 1, -4*xsi, 4*xsi, 4*(1 - xsi - 2*eta)];
    gradN = invB' * gradN_hat;
    % Contribution à la matrice de rigidité
    Kel = Kel + poids * (gradN' * gradN) * detB ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%24
