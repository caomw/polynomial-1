[a, b] = meshgrid(1:9, 1:9);
a = a(:);
b = b(:);
c = zeros(81, 1)
a = [a, b, c];

f = fopen('mesh_9x9.pts', 'wt+');

for i = 1 : size(a, 1)
    fprintf(f, '%2.4f ', a(i, :));
    fprintf(f, '\n');
end
fclose(f);

f = fopen('mesh_9x9.tri', 'wt+');
a = 0 : 80;
a = reshape(a, 9, 9)';
for i = 1 : 8
    for j = 1 : 8
        fprintf(f, '%d ', a(i, j), a(i, j+1), a(i+1, j));
        fprintf(f, '\n');
        fprintf(f, '%d ', a(i, j+1), a(i+1, j+1), a(i+1, j));
        fprintf(f, '\n');
    end
end
fclose(f);


f = fopen('edge_9x9.txt', 'wt+');
for i = 1 : size(E)
   fprintf(f, '%d %d %2.6f ', E(i,1), E(i,2), E(i,3));
   fprintf(f, '\n');
end

