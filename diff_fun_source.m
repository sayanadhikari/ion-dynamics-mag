
function dy = diff_fun_source(~,y)
global gamma_x gamma_z beta

dy=zeros(6,1);


dy(1) = y(2);
dy(2) = exp(y(1)) - y(3);
dy(3) =  ((y(3)/(y(6)^2))*y(2))+(y(3)*gamma_x*(y(5)/(y(6)^2)))+(2*(beta/y(6))*exp(y(1))); 
dy(4) = (gamma_z*(y(5)/y(6)))-(beta*(exp(y(1))/y(3))*(y(4)/y(6)));
dy(5) = (gamma_x)-(gamma_z*(y(4)/y(6)))-(beta*(exp(y(1))/y(3))*(y(5)/y(6)));
dy(6) = -(y(2)/y(6))-(gamma_x*(y(5)/y(6)))-(beta*(exp(y(1))/y(3)));
end