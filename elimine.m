% function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu)
%     % Initialisation
%     tilde_AA = AA; % Copie de la matrice globale
%     tilde_LL = LL; % Copie du vecteur second membre
% 
%     % Parcours des noeuds pour appliquer la condition de Dirichlet
%     for I = 1:length(Refneu)
%         if Refneu(I) == 1 % Si le noeud est sur le bord
%             % Mise à zéro de la ligne et colonne correspondantes
%             tilde_AA(I, :) = 0;
%             tilde_AA(:, I) = 0;
%             % Placement de 1 sur la diagonale
%             tilde_AA(I, I) = 1;
%             % Modification du vecteur second membre
%             tilde_LL(I) = 0;
%         end
%     end
% end
% 

function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu)
    % elimine : Applique les conditions de Dirichlet non homogènes au système linéaire.
    % Entrées :
    %   AA : Matrice globale initiale.
    %   LL : Vecteur membre de droite initial.
    %   Refneu : Vecteur d'indicateurs des nœuds (1 pour les bords Dirichlet).
    %   Coorneu : Matrice des coordonnées des nœuds [x, y].
    % Sorties :
    %   tilde_AA : Matrice globale modifiée.
    %   tilde_LL : Vecteur membre de droite modifié.

    % Trouver les indices des nœuds Dirichlet
    Idx = find(Refneu == 1);

    % Initialisation des matrices
    tilde_AA = AA;
    tilde_LL = LL;

    % Appeler la fonction g pour calculer les valeurs aux nœuds Dirichlet
    g_values = g(Coorneu, Idx);

    % Appliquer les conditions de Dirichlet
    tilde_AA(Idx, :) = 0;          % Mettre les lignes correspondantes à zéro
    tilde_LL(Idx) = g_values;      % Mettre à jour le vecteur membre de droite

    % Fixer les éléments diagonaux à 1
    for i = 1:length(Idx)
        id = Idx(i);
        tilde_AA(id, id) = 1;
    end
end
