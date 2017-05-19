function [Data]=img2patch(I,b)
    r=(b-1)/2;
    [M,N]=size(I);
    Data=zeros(b*b,M*N);
    Data_sub=im2col(I,[b,b],'sliding');
    [x,y]=meshgrid(r+1:N-r,r+1:M-r);
    Idx=sub2ind([M,N],y,x);
    Data(:,Idx)=Data_sub;
    %------------------up boundary---------------------%
    for i=1:r
        for j=1:N
        idx=sub2ind([M,N],i,j);
        i_start=max(1,i-r);
        i_end=min(i_start+b-1,M);
        j_start=max(1,j-r);
        j_end=min(j_start+b-1,N);
        if (j_end-j_start)<b-1
            j_start=j_start-(b-(j_end-j_start))+1;
        end
        Patch=I(i_start:i_end,j_start:j_end);
        Data(:,idx)=Patch(:);
        end
    end
    %------------------bottomp boundary---------------------%
    for i=M-r:M   
        for j=1:N
        idx=sub2ind([M,N],i,j);
        i_start=max(1,i-r);
        i_end=min(i_start+b-1,M);
        j_start=max(1,j-r);
        j_end=min(j_start+b-1,N);
        if (j_end-j_start)<b-1
            j_start=j_start-(b-(j_end-j_start))+1;
        end
        if (i_end-i_start)<b-1
            i_start=i_start-(b-(i_end-i_start))+1;
        end
 
        Patch=I(i_start:i_end,j_start:j_end);
        Data(:,idx)=Patch(:);
        end
    end
     %------------------left boundary---------------------%
    for i=1:M
        for j=1:r
        idx=sub2ind([M,N],i,j);
        i_start=max(1,i-r);
        i_end=min(i_start+b-1,M);
        j_start=max(1,j-r);
        j_end=min(j_start+b-1,N);
        if (j_end-j_start)<b-1
            j_start=j_start-(b-(j_end-j_start))+1;
        end
        if (i_end-i_start)<b-1
            i_start=i_start-(b-(i_end-i_start))+1;
        end
 
        Patch=I(i_start:i_end,j_start:j_end);
        Data(:,idx)=Patch(:);
        end
    end
    %------------------right boundary---------------------%
    for i=1:M
        for j=N-r:N
        idx=sub2ind([M,N],i,j);
        i_start=max(1,i-r);
        i_end=min(i_start+b-1,M);
        j_start=max(1,j-r);
        j_end=min(j_start+b-1,N);
        if (j_end-j_start)<b-1
            j_start=j_start-(b-(j_end-j_start))+1;
        end
        if (i_end-i_start)<b-1
            i_start=i_start-(b-(i_end-i_start))+1;
        end
 
        Patch=I(i_start:i_end,j_start:j_end);
        Data(:,idx)=Patch(:);
        end
    end       
end