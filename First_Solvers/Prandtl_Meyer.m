function angle = Prandtl_Meyer(M)
global gamma
   angle = sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)*(M^2-1)/(gamma+1)))-atan(sqrt(M^2-1));
end