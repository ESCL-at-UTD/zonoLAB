clear
Z = zono([eye(2),ones(2)],ones(2,1));
Z2 = zono([eye(2),2*eye(2)],zeros(2,1));
Z3 = union(Z,Z2);

temp = memZono(Z,'test');
temp2 = memZono(Z2,{'test_2','test_3','test_4','test_5'});
temp3 = memZono(Z3,'test');

temp.dimKeys = 'test';
temp2.dimKeys = {'test_2','test_3'};

% minSum(temp,temp2)
% temp+temp2

% temp3 = and(temp,temp2)

% plot(temp3)


