%Hulei
%20190326
% PLS (NIPALS algorithm) 
function [t,p,q,b,w,r] = pls_NIPALS(X0,Y0, a)
    % a is the number of PLS factors
    [n,m] = size(X0);
    X = zeros(n,m,a);
    Y = zeros(n,size(Y0, 2),a);
    X(:,:,1) = X0;
    Y(:,:,1) = Y0;
    for i = 1:a
        u(:,i) = Y(:,1,i);
        itererr = 1;
        while norm(itererr)>0.0000000001
            w(:,i)=X(:,:,i)'*u(:,i)/norm(X(:,:,i)'*u(:,i));
            w(:,i)= w(:,i)/norm( w(:,i));
            t(:,i)=X(:,:,i)*w(:,i);
            q(:,i)=Y(:,:,i)'*t(:,i)/(t(:,i)'*t(:,i));
            q(:,i)= q(:,i)/norm(q(:,i));
            u(:,i)=Y(:,:,i)*q(:,i);
            itererr=t(:,i)-X(:,:,i)*X(:,:,i)'*u(:,i)/norm(X(:,:,i)'*u(:,i));
        end
        b(:,i) = u(:,i)'*t(:,i)/(t(:,i)'*t(:,i));
        %得分t正交
        p(:,i)=X(:,:,i)'*t(:,i)/(t(:,i)'*t(:,i));
%         t(:,i)=t(:,i)*norm(p(:,i));
%         p(:,i)=p(:,i)/norm(p(:,i));
        X(:,:,i+1)=X(:,:,i)-t(:,i)*p(:,i)';
        Y(:,:,i+1)=Y(:,:,i)-b(:,i)*t(:,i)*q(:,i)'; % y第一层是原始未提取潜变量的y；第二层是提取了第一个潜变量后的残差y，依次类推
    end
    b = diag(b);
%     r=w*inv(p'*w); 
    r =0
end