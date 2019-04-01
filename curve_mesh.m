clear;clc;close all
addpath('../../Project-1/code')

gridname = '../../grid/bump0';
[node,E2N,bdy] = readMesh(gridname);
[I2E,B2E,In,Bn,Area] = processMesh(node,E2N,bdy);

nnode = size(node,1);
E2N_orig = E2N;
bump = @(x) 0.0625*exp(-25*x.^2);
q = 2;
node_num = nnode;
E2N_new = [];
deleted_elem = [];
for ib = 1:size(B2E,1)
    if B2E(ib,3) == 4 % bottom wall
        elem = B2E(ib,1);
        face = B2E(ib,2);
        deleted_elem = [deleted_elem,elem];
        elem_E2N = E2N(elem,:);
        elem_node = node(elem_E2N,:);
        if q == 2
            xy(1,:) = (elem_node(1,:)+elem_node(2,:))/2;
            xy(2,:) = (elem_node(1,:)+elem_node(3,:))/2;
            xy(3,:) = (elem_node(2,:)+elem_node(3,:))/2;
            idx = bface2node(face,q);
%             xy(idx,2) = bump(xy(idx,1));
            node = [node;xy];
            new_elem_E2N = [elem_E2N(1),node_num+1,elem_E2N(2),node_num+2,node_num+3,elem_E2N(3)];
            node_num = node_num + 3;
            E2N_new = [E2N_new;new_elem_E2N];
        elseif q == 3
            xy(1,:) = (elem_node(2,:)-elem_node(1,:))/3 + elem_node(1,:);
            xy(2,:) = 2*(elem_node(2,:)-elem_node(1,:))/3 + elem_node(1,:);
            xy(3,:) = (elem_node(3,:)-elem_node(1,:))/3 + elem_node(1,:);
            xy(6,:) = 2*(elem_node(3,:)-elem_node(1,:))/3 + elem_node(1,:);
            xy(5,:) = (elem_node(3,:) - elem_node(2,:))/3 + elem_node(2,:);
            xy(7,:) = 2*(elem_node(3,:) - elem_node(2,:))/3 + elem_node(2,:);
            xy(4,:) = (elem_node(1,:)+elem_node(2,:)+elem_node(3,:))/3;
            node = [node;xy];
            new_elem_E2N = [elem_E2N(1),node_num+1,node_num+2,elem_E2N(2),...
                node_num+3,node_num+4,node_num+5,...
                node_num+6,node_num+7,elem_E2N(3)];
            node_num = node_num + 7;
            E2N_new = [E2N_new;new_elem_E2N];
        end
    end
end
% E2N(deleted_elem,:) = [];
writeCurvedMesh(strcat(gridname,'_',num2str(q),'_linear'),q,node,E2N_orig,E2N,E2N_new,deleted_elem);
% plot(node(:,1),node(:,2),'.')
% axis equal
% plotgri('../../grid/bump0_3');

function [xy,idx2] = getNode(node,E2N,elem,face)
    idx = E2N(elem,:);
    xy = zeros(2,2);
    j = 1;
    for i = 1:3
        if i ~= face
            xy(j,:) = node(idx(i),:);
            idx2(j) = i;
            j = j + 1;
        end
    end
end

function [y] = bface2node(bface,q)
    if q == 2
        if bface == 1
            y = 3;
        elseif bface == 2
            y = 2;
        elseif bface == 3
            y = 1;
        end
    elseif q == 3
        if bface == 1
            y = [5,7];
        elseif bface == 2
            y = [3,6];
        elseif bface == 3
            y = [1,2];
        end
    end
end