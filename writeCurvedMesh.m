function writeCurvedMesh(name,q,node,E2N_orig,E2N,E2N2,qlist)
    nNode = size(node,1);
    nElemTot = size(E2N_orig,1);
    nElem = size(E2N,1);
    nElem2 = size(E2N2,1);
    nBGroup = 4;
    xmin = -1.5;
    xmax = 1.5;
    ymax = 0.8;
    tol = 1e-8;
    
    fid = fopen(sprintf('%s.gri',name), 'w');
    
    fprintf(fid, '%8d %8d %8d\n', nNode, nElemTot, 2);
    for i = 1:nNode
        x = node(i,1);
        y = node(i,2);
        fprintf(fid, '%25.16E %25.16E\n', x, y);
    end
    fprintf(fid, '%d\n', nBGroup); % # boundary face groups
    
    
    [~,BE] = edgeHash(E2N_orig);
    nbedge = size(BE,1);
    top = zeros(nbedge); % We just initialize the maximum size needed
    left = zeros(nbedge);
    right = zeros(nbedge);
    bottom = zeros(nbedge);
    nodeL = zeros(nbedge,2);
    nodeR = zeros(nbedge,2);
    nodeT = zeros(nbedge,2);
    nodeB = zeros(nbedge,2);
    nL = 0;
    nR = 0;
    nT = 0;
    nB = 0;
    for i = 1:nbedge
        n1 = BE(i,1);
        n2 = BE(i,2);
        if abs(node(n1,1)-xmin) < tol && abs(node(n2,1)-xmin) < tol
            nL = nL + 1;
            left(nL) = i;
            nodeL(nL,:) = [n1,n2];
        elseif abs(node(n1,1)-xmax) < tol && abs(node(n2,1)-xmax) < tol
            nR = nR + 1;
            right(nR) = i;
            nodeR(nR,:) = [n1,n2];
        elseif abs(node(n1,2)-ymax) < tol && abs(node(n2,2)-ymax) < tol
            nT = nT + 1;
            top(nT) = i;
            nodeT(nT,:) = [n1,n2];
        else
            nB = nB + 1;
            bottom(nB) = i;
            nodeB(nB,:) = [n1,n2];
        end
    end
    
    fprintf(fid, '%d %8d Left\n', nL, 2);
    for i = 1:nL
        fprintf(fid, '%8d %8d\n', nodeL(i,1),nodeL(i,2));
    end
    
    fprintf(fid, '%d %8d Right\n', nR, 2);
    for i = 1:nR
        fprintf(fid, '%8d %8d\n', nodeR(i,1),nodeR(i,2));
    end
    
    fprintf(fid, '%d %8d Top\n', nT, 2);
    for i = 1:nT
        fprintf(fid, '%8d %8d\n', nodeT(i,1),nodeT(i,2));
    end
    
    fprintf(fid, '%d %8d Bottom\n', nB, 2);
    for i = 1:nB
        fprintf(fid, '%8d %8d\n', nodeB(i,1),nodeB(i,2));
    end
    
    fprintf(fid, '%d 1 TriLagrange\n', nElem);
    for i = 1:nElem
        fprintf(fid, '%8d %8d %8d\n', E2N(i,1),E2N(i,2),E2N(i,3));
    end
    
    fprintf(fid, '%d %d TriLagrange\n', nElem2, q);
    for i = 1:nElem2
        for iq = 1:(q+1)*(q+2)/2
            fprintf(fid, '%8d', E2N2(i,iq));
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '%d qlist\n', numel(qlist));
    for i = 1:numel(qlist)
        fprintf(fid, '%8d\n', qlist(i));
    end
    
    fclose(fid);
%end