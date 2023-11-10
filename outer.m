function [C]=outer(A,B)
sizeA = size(A);
sizeB = size(B);
C     = zeros([sizeA, sizeB]);
if ndims(A) == 2 && ndims(B) == 2
    for b2 = 1:sizeB(2)
        for b1 = 1:sizeB(1)
            for a2 = 1:sizeA(2)
                for a1 = 1:sizeA(1)
                    C(a1, a2, b1, b2) = A(a1, a2) * B(b1, b2);
                end
            end
        end
    end
    C=squeeze(C);
elseif ndims(A) == 2 && ndims(B) == 3
    for b3=1:sizeB(3)
        for b2 = 1:sizeB(2)
            for b1 = 1:sizeB(1)
                for a2=1:sizeA(2)
                    for a1 = 1:sizeA(1)
                        C(a1, a2, b1, b2, b3) = A(a1) * B(b1, b2,b3);
                    end
                end 
            end
        end
    end
    C=squeeze(C);
else
    error('Error.');
end
end
