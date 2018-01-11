function [ Val ] = eval_at_point( Ord,x,uh,x0 )

  xi  = (x0-x(ne))/(x(ne+1)-x(ne))*2+(-1);
  un  = basis_1d(Ord,xi);
  Val = uh(ne,:)*un;

end
