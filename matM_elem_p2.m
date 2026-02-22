function [Mel] = matM_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de masse elementaire en P2 Lagrange
%
% SYNOPSIS [Mel] = matM_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de masse
% -------------------------------
% Initialisation

Mel =zeros(6, 6);
B = [x2 - x1, x3 - x1;
    y2 - y1, y3 - y1];
if (abs(det(B)) <= eps)
    error('l aire d un triangle est nulle!!!');
end;
detB = abs(det(B));

% Points et poids de quadrature
S_hat = [0.0915762135098, 0.0915762135098;
    0.8168475729805, 0.0915762135098;
    0.0915762135098, 0.8168475729805;
    0.1081030181681, 0.4459484909160;
    0.4459484909160, 0.1081030181681;
    0.4459484909160, 0.4459484909160];
poids = [0.05497587183, 0.05497587183, 0.05497587183, 0.1116907948, 0.1116907948, 0.1116907948];

% Boucle sur les points de quadrature

for q = 1:length(poids)
    % Coordonnées des points de quadrature
    lambda1 = 1 - S_hat(q, 1) - S_hat(q, 2); % Coordonnée barycentrique 1
    lambda2 = S_hat(q, 1);                   % Coordonnée barycentrique 2
    lambda3 = S_hat(q, 2);                   % Coordonnée barycentrique 3
    % On exprime w sous forme d'un vecteur afin d'alléger le calcul
    w = [lambda1 * (2 * lambda1 - 1);
        lambda2 * (2 * lambda2 - 1);
        lambda3 * (2 * lambda3 - 1);
        4 * lambda1 * lambda2;
        4 * lambda2 * lambda3;
        4 * lambda1 * lambda3];
    Mel = Mel + poids(q) * (w * w') * detB ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%