%构建RMSE函数
function result=RMSE(A,B)
   M=length(A);
   a=0;sum=0;
   for i=1:M
      for j=1:M
          a=(A(i,j)-B(i,j)).^2;
          sum=sum+a;
      end
   end
      result=sqrt(sum/(M*M));%注意由于本次的M=N=512
end



