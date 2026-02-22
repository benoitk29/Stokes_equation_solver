% =====================================================
%
% une routine pour la mise en oeuvre des EF P2 Lagrange
% pour l'equation de Laplace suivante, avec conditions de Neumann
%
% | -\Delta u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================
close all;
clear all;
h_valeurs = [0.2, 0.1, 0.05, 0.025];  % Les différents h
% fichiers_mesh = {'geomRectangle_h02.msh', 'geomRectangle_h01.msh', 'geomRectangle_h005.msh', 'geomRectangle_h0025.msh'};
fichiers_mesh = {'geomRectangle_h02.msh','geomRectangle_h01.msh', 'geomRectangle_h005.msh','geomRectangle_h0025.msh'};

% Initialisation des erreurs
L2_erreurs = zeros(length(h_valeurs), 1);
H1_erreurs = zeros(length(h_valeurs), 1);
L2_exact =  zeros(length(h_valeurs), 1);
L2_relative =zeros(length(h_valeurs), 1);
H1_exact = zeros(length(h_valeurs), 1);
H1_relative=zeros(length(h_valeurs), 1);

% lecture du maillage et affichage
% ---------------------------------
for idx = 1:length(h_valeurs)
    nom_maillage = fichiers_mesh{idx};
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]= lecture_msh_ordre2(nom_maillage);
    affichemaillage_ordre2(nom_maillage);
    xlabel('x');
    ylabel('y');

    % ----------------------
    % calcul des matrices EF
    % ----------------------

    % declarations
    % ------------
    KK = sparse(Nbpt,Nbpt); % matrice de rigidite
    MM = sparse(Nbpt,Nbpt); % matrice de masse
    LL = zeros(Nbpt,1);     % vecteur second membre

    % boucle sur les triangles
    % ------------------------
    for l=1:Nbtri

        % Coordonnees des sommets du triangles
        %De même que dans le TP1
        S1 = Coorneu(Numtri(l, 1), :);
        S2 = Coorneu(Numtri(l, 2), :);
        S3 = Coorneu(Numtri(l, 3), :);

        % Calcul des matrices elementaires du triangle l
        Kel=matK_elem_p2(S1, S2, S3);
        Mel=matM_elem_p2(S1, S2, S3);

        % On fait l'assemblage des matrices globales
        % De même que dans le TP1
        for i=1:6
            I= Numtri(l,i);
            for j=1:6
                J=Numtri(l,j);
                MM(I,J)= MM(I,J)+Mel(i,j);
                KK(I,J)= KK(I,J)+Kel(i,j);
            end %boucle j
        end %boucle i

    end % boucle l

    AA=MM+KK;
    % Calcul du second membre
    % -------------------------

    % utiliser la routine f.m
    FF = zeros(Nbpt,1);
    FF = f(Coorneu(:,1),Coorneu(:,2));
    LL = MM*FF;


    % Resolution du systeme lineaire
    % ----------

    UU = AA\LL;

    % visualisation
    % -------------
    affiche_ordre2(UU, Numtri, Coorneu,sprintf('Neumann - %s', nom_maillage));
    xlabel('x');
    ylabel('y');

    validation = 'oui';
    % validation
    % ----------
    if strcmp(validation,'oui')
        %Affichage de la solution exacte
        UU_exact = 3 * cos(pi * Coorneu(:, 1)) .* cos(2 * pi * Coorneu(:, 2));
        affiche_ordre2(UU_exact, Numtri, Coorneu, sprintf('Neumann solution exacte- %s', nom_maillage));
        xlabel('x');
        ylabel('y');

        % Calcul de l'erreur L2
        L2_erreurs(idx) = sqrt((UU - UU_exact)' * MM * (UU - UU_exact));
        L2_exact(idx) = sqrt(UU_exact' * MM * UU_exact);
        L2_relative(idx) = L2_erreurs(idx)./L2_exact(idx);


        % Calcul de l'erreur H1
        H1_erreurs(idx) = sqrt((UU - UU_exact)' * KK * (UU - UU_exact));
        H1_exact(idx) = sqrt(UU_exact' * KK * UU_exact);
        H1_relative(idx) = H1_erreurs(idx)./H1_exact(idx);
        % attention de bien changer le terme source (dans FF)
    end
end

%Affichage optionnel des valeurs de h, erreurs L2/H1 relatives et absolues
%dans la fenêtre de commande de Matlab
fprintf('Pas de h: %f\n', h_valeurs);
fprintf('Erreur L2: %f\n', L2_erreurs);
fprintf('Erreur relative L2: %f\n', L2_relative);
fprintf('Erreur H1: %f\n', H1_erreurs);
fprintf('Erreur relative H1: %f\n', H1_relative);
log_h = log10(1 ./ h_valeurs);

% Régression linéaire pour L2 et H1
poly_L2 = polyfit(log_h, log10(L2_relative), 1);
poly_H1 = polyfit(log_h, log10(H1_relative), 1);
regression_L2 = polyval(poly_L2, log_h);
regression_H1 = polyval(poly_H1, log_h);

%---------------------------Partie graphique---------------------------
figure;
% Tracé pour l'erreur L2 et H1 dans le même graphique
plot(log_h, log10(L2_relative), '-o', 'DisplayName', 'norme L2', 'Color', 'r', 'Marker', '*');
hold on;
plot(log_h, log10(H1_relative), '-o', 'DisplayName', 'seminorme H1', 'Color', 'b', 'Marker', 'o');
plot(log_h, regression_L2, '--', 'DisplayName', 'Courbe de regression de L2');
plot(log_h, regression_H1, '--', 'DisplayName', 'Courbe de regression de H1');

legend('show', 'Location', 'best');
xlabel('log(1/h)');
ylabel('log(Erreur relative)');
title('Erreur L2 et H1 en fonction de h');
legend('show');
hold off;

% Affichage des pentes de régression (ordre de convergence)
fprintf('Pente de la régression pour l''erreur L2 : %f\n', poly_L2(1));
fprintf('Pente de la régression pour l''erreur H1 : %f\n', poly_H1(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%24
