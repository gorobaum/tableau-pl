function [ind, x] = simplex(A,b,c,m,n,print)
  [ind, x] = preparesimplex(A,b,c,m,n,print)
endfunction

function [ind, x] = preparesimplex(A,b,c,m,n,print)
  auxb = [];
  for i = 1:m
    if ( b(i) < 0.0 )
      auxb = [auxb;-b(i)];
      auxA(i,:) = -A(i,:);
    else
      auxb = [auxb;b(i)];
      auxA(i,:) = A(i,:);
    endif
  endfor
  cf1 = [] ;
  for i = 1:n
    cf1 = [cf1;0.0];
  endfor
  for j = (i+1):m+n
    indbase(j-i) = j;
    cf1 = [cf1;1.0];
  endfor
  
  I = eye(m);
  auxA = [auxA,I];
  auxc(1) = -(base(cf1,indbase))'*auxb;
  auxc = [auxc(1),cf1'-(base(cf1,indbase))'*auxA];
  auxA = [auxb,auxA];
  auxA = [auxc;auxA];
  
  disp("Fase 1")
  [ind, x, B, indb, m] = runsimplex(auxA,auxb,auxc,indbase,m,m+n,print)
  B = B(:,1:n+1);
  [B, indb, m] = removeslackformbase(B,m,n,indb)
  if (ind == -1 )
    ind = 1;
  else
    disp("Fase 2")
    auxc(1) = -(base(c,indb))'*B(2:m+1,1)
    auxc = [auxc(1),c'-(base(c,indb))'*B(2:m+1,2:n+1)]
    B(1,:) = auxc
    disp("Phase 2 hativaaa")
    [ind, x, B] = runsimplex(B,b,c,indb,m,n,print)
  endif
endfunction

function [ind, x, B, indb, m] = runsimplex(A,b,c,indbase,inm,n,print)
  
  stop = false
  while (!stop)
    disp("Novo passo")
    A;
    j = 2;
    m = inm;
    while ( j <= length(A(1,:)) && A(1,j) >= 0.0 )
      j++;
    endwhile
    if ( j == length(A(1,:))+1 )
      stop = true
      ind = 0
    else
      j
      u = A(:,j)
      ratio = inf;
      for i = 2:length(u)
        if ( u(i) > 0.0 && ratio > A(i,1)/u(i) )
          l = i;
          ratio = A(i,1)/u(i);
        endif
      endfor
      if ( i == length(u)+1 && ratio == inf)
        stop = true
        ind = -1
      else
        indbase(l-1) = j-1
        for i = 1:length(u)
          if ( i != l && u(i) != 0.0 )
            A(i,:) -= (u(i)/u(l))*A(l,:);
          endif
        endfor
        A(l,:) /= u(l);
      endif
    endif
  endwhile
  disp("Final do simplex")
  indb = indbase;
  x = A(2:m+1,1);
  B = A;
endfunction

function [B, ind, m] = removeslackformbase(A,inm,n,indbase)
  disp("Retirada da base")
  B = A
  ind = indbase
  length(indbase)
  m = inm
  for i = 2:m
    for j = m+1:m+n
      if (indbase(i-1) == j )
        k = 1;
        while ( k <= length(A(i,:)) && A(i,k) == 0.0 )
          k++;
        endwhile
        if ( k == length(A(i+1,:))+1 )
          i
          B = [B(1:i-1,:);B(i+1:m+1,:)];
          ind = [ind(1:i-2),ind(i:m)];
          m--;
        else
        endif
      endif
    endfor
  endfor
  disp("Fim da retirada")
endfunction

function ret = base(vec, b)
  j = 1;
  ret = [];
  for i = 1:length(vec)
    if ( j <= length(b) && i == b(j) )
      ret = [ret;vec(i)];
      j++;
    endif
  endfor
endfunction
