clear
Z = zono(eye(2),ones(2,1));

temp = memZono(Z,'test');
temp2 = memZono(2*Z,{'test_2','test_3'});

% minSum(temp,temp2)
% temp+temp2