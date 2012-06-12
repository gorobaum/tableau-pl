function [ind, x] = simplex(A,b,c,m,n,print)
  [B, r] = passo1(A,b,c,m,n,print)
endfunction

function [B, r] = passo1(A,b,c,m,n,print)
  for i = 1:length (b)
    if ( b(i) < 0)
      b(i) = -1 * b(i)
    endif
  endfor
  B = eye(m)
  for i = 1:n
    x(i) = 0
  endfor
  [ind, x] = realsimplex(A,B,b,c,m,n,print)
endfunction

function [ind, x] =  realsimplex(A,B,b,c,m,n,print)
  binv = inv(B)
  j = 0
  
  while ( c(j) > 0.0 )
    j++
  endwhile
  if ( j == length(c) )
    ind = 0
    stop = true
  endif
  
  u = binv*A(:,j)
  i = 0
  while ( u(i) < 0.0 )
    i++
  endwhile
  if ( i == length(u) )
    ind = -1
    stop = true
  endif
endfunction
