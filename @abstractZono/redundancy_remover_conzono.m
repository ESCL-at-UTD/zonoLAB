function [G,c,A,b,eff_r]= redundancy_remover_conzono(G,c,A,b)


zeroCols=find(all([G;A]==0,1));
G(:,zeroCols)=[];
A(:,zeroCols)=[];
zerorows=find(all([A,b]==0,2));
A(zerorows,:)=[];
b(zerorows)=[];

E=bound_e_LP_calculator(G,A,b);

equal_indices = find(E(:, 1) == E(:, 2));
c = c + sum(G(:, equal_indices) .*E(equal_indices, 1)', 2);
G(:,equal_indices)=[];
b=b-sum(A(:,equal_indices).*E(equal_indices,1)',2);
A(:,equal_indices)=[];

R=bound_r_LP_calculator(G,A,b);

% [E,R]=refine_bounds_vectorized_function(A,b);
% 
% equal_indices = find(E(:, 1) == E(:, 2));
% 
% c = c + sum(G(:, equal_indices) .* E(equal_indices, 1)', 2);
% G(:,equal_indices)=[];
% b=b-sum(A(:,equal_indices).*E(equal_indices,1)',2);
% A(:,equal_indices)=[];
% %obj=conZono(G,c,A,b);
% % R(equal_indices,:)=[];
% % E(equal_indices,:)=[];
% 
% [E,R]=refine_bounds_vectorized_function(A,b);
% 
% E_m=1/2*(E(:,1)+E(:,2));
% E_r=1/2*(E(:,2)-E(:,1));     
% G_new=G*diag(E_r);
% c_new=c+G*E_m;
% A_new=A*diag(E_r);
% b_new=b-A*E_m;
% 
% G=G_new;
% c=c_new;
% A=A_new;
% b=b_new;
% 
% R(:,1)=(R(:,1)-E_m)./E_r;
% R(:,2)=(R(:,2)-E_m)./E_r;

eff_r=find((R(:,1) >= -1) & (R(:,2) <= 1));
% close_indices=abs(E(eff_r,1)-E(eff_r,2))<1e-2;
% eff_r(close_indices)=[];

end
