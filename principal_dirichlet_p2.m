% Valeurs de h et noms des fichiers de maillage correspondants
h_values = [0.2, 0.1, 0.05, 0.025];  % Les différents h
mesh_files = {'geomRectangle_h02.msh', 'geomRectangle_h01.msh', 'geomRectangle_h005.msh',  'geomRectangle_h0025.msh'};

% Initialisation des erreurs
L2_errors = zeros(length(h_values), 1);
H1_errors = zeros(length(h_values), 1);

for idx = 1:length(h_values)
    nom_maillage = mesh_files{idx};
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]= lecture_msh_ordre2(nom_maillage);

    % Assemblage des matrices
    MM = sparse(Nbpt, Nbpt); % Matrice de masse
    KK = sparse(Nbpt, Nbpt); % Matrice de rigidité
    LL = zeros(Nbpt, 1);     % Vecteur second membre

    for l = 1:Nbtri
        S1 = Coorneu(Numtri(l, 1), :);
        S2 = Coorneu(Numtri(l, 2), :);
        S3 = Coorneu(Numtri(l, 3), :);

        Kel = matK_elem_p2(S1, S2, S3);
        Mel = matM_elem_p2(S1, S2, S3);

        for i = 1:6
            I = Numtri(l, i);
            for j = 1:6
                J = Numtri(l, j);
                MM(I, J) = MM(I, J) + Mel(i, j);
                KK(I, J) = KK(I, J) + Kel(i, j);
            end
        end
    end

    AA = MM + KK; % Matrice du système
    FF = f_dirichlet(Coorneu(:, 1), Coorneu(:, 2));
    LL = MM * FF; % Second membre sans condition de bord

    % Pseudo-élimination
    [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu); 

    % Résolution du système linéaire
    UU = tilde_AA \ tilde_LL;

    % visualisation
    % -------------
    affiche_ordre2(UU, Numtri, Coorneu,sprintf('Dirichlet - %s', nom_maillage));

    validation = 'oui';
    % validation
    % ----------
    if strcmp(validation,'oui')
        UU_exact = 3*sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
        % A COMPLETER
        L2_errors(idx) = sqrt((UU - UU_exact)' * MM * (UU - UU_exact));
        L2_exact(idx) = sqrt(UU_exact' * MM * UU_exact);
        L2_relative(idx) = L2_errors(idx)./L2_exact(idx);
        fprintf('Erreur L2: %f\n', L2_errors);
        fprintf('Erreur relative L2: %f\n', L2_relative);
        % Calcul de l erreur H1
        H1_errors(idx) = sqrt((UU - UU_exact)' * KK * (UU - UU_exact));
        H1_exact(idx) = sqrt(UU_exact' * KK * UU_exact);
        H1_relative(idx) = H1_errors(idx)./H1_exact(idx);
        fprintf('Erreur H1: %f\n', H1_errors);
        fprintf('Erreur relative H1: %f\n', H1_relative);
        % attention de bien changer le terme source (dans FF)
    end
end

log_h = log10(1 ./ h_values);

% Régression linéaire pour L2 et H1
poly_L2 = polyfit(log_h, log10(L2_relative), 1);
poly_H1 = polyfit(log_h, log10(H1_relative), 1);
regression_L2 = polyval(poly_L2, log_h);
regression_H1 = polyval(poly_H1, log_h);

figure;
% Tracé pour l'erreur L2 et H1 dans le même graphique
plot(log_h, log10(L2_relative), '-o', 'DisplayName', 'Norme L2', 'Color', 'r', 'Marker', '*'); % Courbe L2 en rouge
hold on;
plot(log_h, log10(H1_relative), '-o', 'DisplayName', 'Semi Norme H1', 'Color', 'b', 'Marker', 'o'); % Courbe H1 en bleu
plot(log_h, regression_L2, '--', 'DisplayName', 'Regression Norme L2');
plot(log_h, regression_H1, '--', 'DisplayName', 'Regression Norme H1');

legend('show', 'Location', 'best');
xlabel('log(1/h)');
ylabel('log(Erreur)');
title('Erreur L2 et H1 en fonction de h');
legend('show');
hold off;

% Affichage des pentes de régression (ordre de convergence)
fprintf('Pente de la régression pour l''erreur L2 : %f\n', poly_L2(1));
fprintf('Pente de la régression pour l''erreur H1 : %f\n', poly_H1(1));