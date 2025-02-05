function [G,c,A,b,eff_r]= redundancy_remover_conzono(G,c,A,b)

zeroCols=find(all([G;A]==0,1));
G(:,zeroCols)=[];
A(:,zeroCols)=[];
zerorows=find(all([A,b]==0,2));
A(zerorows,:)=[];
b(zerorows)=[];

[E,R]=refine_bounds_function(A,b);

equal_indices = find(E(:, 1) == E(:, 2));

c = c + sum(G(:, equal_indices) .* E(equal_indices, 1)', 2);
G(:,equal_indices)=[];
b=b-sum(A(:,equal_indices).*E(equal_indices,1)',2);
A(:,equal_indices)=[];
%obj=conZono(G,c,A,b);
R(equal_indices,:)=[];

eff_r=find((R(:,1) >= -1) & (R(:,2) <= 1));

end
