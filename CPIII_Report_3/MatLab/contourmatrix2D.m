function contourmatrix2D(A,b,x1,x2,y1,y2)
    X=x1:(x2-x1)/300:x2;
    Y=y1:(y2-y1)/300:y2;
    Z=zeros(length(X),length(Y));
    r=@(x,y)norm(A*[x,y]'-b);
    for ii=1:length(X)
        for jj=1:length(Y)
            Z(ii,jj)=r(X(ii),Y(jj));
        end
    end
    figure;
    contourf(X,Y,Z);
    set(gca,'Fontsize',20);
    % saveas(gcf,'mat1.eps','epsc');
end