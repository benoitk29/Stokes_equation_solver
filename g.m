function g_values = g(Coorneu, Idx)
    % Définir l'expression de g(x, y)
    % g_expression = @(x, y) 3 * cos(pi * x) .* cos(2 * pi * y);
    g_expression = @(x,y) (0 * x) .* (0 * y);

    % Extraire les coordonnées des nœuds spécifiés
    x = Coorneu(Idx, 1); % Coordonnées x
    y = Coorneu(Idx, 2); % Coordonnées y

    % Calculer les valeurs de g
    g_values = g_expression(x, y);
end
