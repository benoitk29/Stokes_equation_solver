function [Eel] = matE_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice elementaire du bloc rectangulaire (p, dv1/dx)
%
% SYNOPSIS [Eel] = matE_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Eel matrice elementaire rectangulaire (matrice 6x3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Points et poids de quadrature
S_hat = [1/6, 1/6;
    2/3, 1/6;
    1/6, 2/3];
poids = 1/6;

% Calcul de la matrice elementaire du bloc rectangulaire (p, dv1/dx)
% ------------------------------------------------------------------
Eel = zeros(6,3);
%Traitement de Bl pour calculer la transposée de son inverse utilisée dans
% le calcul de quadrature, tout comme son déterminant
B = [x2 - x1, x3 - x1;
    y2 - y1, y3 - y1];
if (abs(det(B)) <= eps)
    error('l aire d un triangle est nulle!!!');
end;
detB = abs(det(B));
invB = inv(B);
invBT=invB';
L1 = invBT(1,:); % Ici, nous prenons seulement la première ligne de l'inverse transposée de Bl


for q = 1:size(S_hat, 1)
    xsi = S_hat(q, 1);
    eta = S_hat(q, 2);
        % On reprends le TP1 pour les fonctions de base P1, on remarque que le
    % D du TP1 est égal au detB. Dans le TP1, les fonctions de base étaient
    % égales aux coordonnées barycentriques.
    %Donc, on peut écrire que
    %     lambda1 = 1 - S_hat(q, 1) - S_hat(q, 2); % Coordonnée barycentrique 1
    %     lambda2 = S_hat(q, 1);                  % Coordonnée barycentrique 2
    %     lambda3 = S_hat(q, 2);                  % Coordonnée barycentrique 3
    % Pour simplifier alléger le code, nous mettons ceci sous forme d'un
    % vecteur w_P1 (dans l'optique de faire un produit scalaire
    w_P1= [1 - S_hat(q, 1) - S_hat(q, 2);S_hat(q, 1);S_hat(q, 2)]; %(surement que que ce n'est y1 mais plutot y1_chapeau, et que x)



    %Après avoir calculé chaque gradient de chaque fonction de bases wi, on
    %fait des 6 gradients trouvés une matrice 6x6 afin de calcul la matrice
    %élémentaire via un produit de matrices
    gradN_hat = [-3 + 4*xsi + 4*eta,            4*xsi - 1,         0,        4*(1 - 2*xsi - eta),     4*eta,        -4*eta;
        -3 + 4*xsi + 4*eta, 0, 4*eta - 1, -4*xsi, 4*xsi, 4*(1 - xsi - 2*eta)];
    wb = w_P1*L1 ;
    % Contribution à la matrice Eel
    Eel = Eel + poids *(-(wb * gradN_hat))'* detB ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%24
