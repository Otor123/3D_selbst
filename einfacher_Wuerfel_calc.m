knoten = [0,0,0;0,1,0;1,0,0;1,1,0;0,0,1;0,1,1;1,0,1;1,1,1];
staebe = [1,2;1,3;1,5;2,4;2,6;3,4;3,7;4,8;5,6;5,7;6,8;7,8;1,8];


% Plotten der aktualisierten Stäbe
figure;
plot3(knoten(:, 1), knoten(:, 2), knoten(:, 3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('X-Achse');
ylabel('Y-Achse');
zlabel('Z-Achse');
title('Knoten');
hold on;
for i = 1:size(staebe, 1)
    k1 = staebe(i, 1);
    k2 = staebe(i, 2);
    x = [knoten(k1, 1), knoten(k2, 1)];
    y = [knoten(k1, 2), knoten(k2, 2)];
    z = [knoten(k1, 3), knoten(k2, 3)];
    plot3(x, y, z, 'r', 'LineWidth', 2);
end
hold off;

% Die Variable "staebe" enthält jetzt die Verbindungen der diagonale Verstrebungen.
disp('Verbindungen (Staebe) mit diagonalen Verstrebungen:');
disp(staebe);