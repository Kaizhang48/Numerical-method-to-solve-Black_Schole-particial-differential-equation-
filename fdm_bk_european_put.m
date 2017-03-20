function V=fdm_bk_european_put(M, N, K, T, r, q, sigma)
% time: N; time step size:k;
% space : M; space step size:h;
% K: strike price
% r interest rate
% q: divident rate
% sigma: volitility
  a=0;
  b=3.*K;
  h=(b-a)./M;
  k=(T-0)./N;
  s=[a,a+(1:M).*h];
  % t=[0,0+(1:N).*k];
  A(1:M)=-(power(sigma,2).*power(s(1:M),2))+(r-q).*h.*s(1:M);
  B(1:M)=2.*power(h,2).*((1/k)+r)+2.*power(sigma,2).*power(s(1:M),2);
  C(1:M-1)=power(sigma,2).*power(s(1:M-1),2)+(r-q).*h.*s(1:M-1);
  Q=diag(B,0)+diag((A(2:M)),-1)+diag((-C),1);
  Q(1,2)= Q(1,2)+A(1);
  V(:,M+1)=0;
  V(1,:)=max(K-s(1:M+1),0);
  for i=2:N+1
    r(1:M,1)=[((2.*power(h,2)/k).*V(i-1,1))-2.*A(1).*h;(2.*power(h,2)/k).*V(i-1,2:M)'];
    V(i,1:M)=(Q\r)';
  end
  % V=V(N+1,:)';
  V=V;
end

