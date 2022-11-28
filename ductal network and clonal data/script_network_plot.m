% ----------------------------------------------------------------------- %
% Cell fate specification in the developing mouse salivary gland is 
% mediated through a cellular hierarchy controlled by Notch and Kras 
% signaling
%
% Chatzeli et al.
% Developmental Cell (2022)
% ----------------------------------------------------------------------- %
%
% This script loads "edge list" and "node position" together with "clone
% position" files to reconstruct the ductal network and clone distribution 
% in the salivary gland.
%
% Inputs:
% edge_list: Array of size number_of_links*4, containing the source node id 
% (column 1), target node id (column 2), and distance from source to target 
% (column 3), and width of the corresponding duct (column 5).
%
% node_positions: Array of size number_of_nodes*4. Column 1: node id. 
% Columns 2-4: coordinates (x,y,z) for each node, respectively.
%
% clone_positions: Array where every row shows the information for a single
% cell. Column 1: number of the clone to which the cell corresponds.  
% Columns 2: type of cell coordinates (x,y,z) for each node, respectively.
% 
% Cell type info:
% clone info: Clone number | index	| x	| y	| z	| cell type
% cell types (for cases where n_cell_types == 3):
% acinar distal     1
% acinar proximal	2
% ductal            3
% cell types (for cases where n_cell_types == 8)
% acinar high Mist1	1
% acinar low Mist1	2
% ductal intercallated	3
% ductal intralobular luminal	4
% ductal intralobular basal	5
% ductal interlobular luminal	6
% ductal inerloboular basal	7
% myoepithelial	8
% ----------------------------------------------------------------------- %
clc; clear all; close all;

% provide path to "data" directory
data_dir = 'E:\Documents\Ignacio\Salivary gland\scripts_cell_fate_salivary_gland_git\data';
remove_ductal_clones = 1; % set to 1 to remove clones with only ductal cells
% select dataset to plot
% Options:
% 'E14.5 control'
% 'E16.5 control'
% 'E18.5 control'
% 'E16.5 kras'
dataset = 'E18.5 control'; 

main_folder_path = {}; file_name_bare = {}; leg = {}; n_cell_types = []; dr = [];

switch dataset
    case 'E14.5 control'
        main_folder_path{end+1} = [data_dir,'\E14_5\control\LCA10_3g 2020.08.02'];n_cell_types(end+1) = 3;
        leg{end+1} = 'E14.5 confetti 1'; 
        dr(end+1,:) = [0.5679, 0.5679, 2]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E14_5\control\LCA10_3g 2020.08.02c'];n_cell_types(end+1) = 3;
        leg{end+1} = 'E14.5 confetti 2'; 
        dr(end+1,:) = [0.5678, 0.5678, 2]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E14_5\control\LCA10_3g 2020.08.03']; n_cell_types(end+1) = 3;
        leg{end+1} = 'E14.5 confetti 3';
        dr(end+1,:) = [0.5679, 0.5679, 2]; % um/pixel

    case 'E16.5 control'
        main_folder_path{end+1} = [data_dir,'\E16_5\LCA152.1e 20200724']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E16.5 confetti 1'; 
        dr(end+1,:) = [0.5678, 0.5678, 2]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E16_5\LCA152.1e 20200729']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E16.5 confetti 2'; 
        dr(end+1,:) = [0.5678, 0.5678, 2]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E16_5\LCA152.1e 20200802']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E16.5 confetti 3'; 
        dr(end+1,:) = [0.5678, 0.5678, 3]; % um/pixel


    case 'E18.5 control'
        main_folder_path{end+1} = [data_dir,'\E18_5\control\E14.5_E18.5\LCA103_1c']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 confetti 1'; 
        dr(end+1,:) = [0.5682, 0.5682, 1]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E18_5\control\E14.5_E18.5\LCA103.1c 2020.03.02']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 confetti 2'; 
        dr(end+1,:) = [0.5678, 0.5678, 2]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E18_5\control\E14.5_E18.5\LCA152.1j 2020.10.09']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 confetti 3';
        dr(end+1,:) = [0.5678, 0.5678, 3]; % um/pixel

    case 'E16.5 kras'
        main_folder_path{end+1} = [data_dir,'\E18_5\r2Kras\LCA98.1h 2020.03.12']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 kras 1'; 
        dr(end+1,:) = [0.5678,0.5678,1.4992]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E18_5\r2Kras\LCA98.1h 2020.03.15']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 kras 2'; 
        dr(end+1,:) = [0.5678,0.5678,2.0014]; % um/pixel
        main_folder_path{end+1} = [data_dir,'\E18_5\r2Kras\LCA98.1h 2020.03.17']; n_cell_types(end+1) = 8;
        leg{end+1} = 'E18.5 kras 3'; 
        dr(end+1,:) = [0.5678,0.5678,2.0014]; % um/pixel
end

%%

% 5. source_node of the network
source_node = 1; 

n_files = length(main_folder_path);

for n_file = 1:n_files
an_folder = [main_folder_path{n_file},'/','out_analysis/'];

edge_list{n_file} = load([an_folder,'edge_list.dat']);
node_positions{n_file} = load([an_folder,'node_positions.dat']);

edge_list_clean{n_file} = load([an_folder,'edge_list_clean.dat']);
node_positions_clean{n_file} = load([an_folder,'node_positions_clean.dat']);

edge_list_large{n_file} = load([an_folder,'edge_list_large.dat']);
node_positions_large{n_file} = load([an_folder,'node_positions_large.dat']);
source_node_large = 1; 

edge_list_large_all_nodes{n_file} = load([an_folder,'edge_list_large_all_nodes.dat']);
node_positions_large_all_nodes{n_file} = load([an_folder,'node_positions_large_all_nodes.dat']);

% clone positions columns: 
if exist([an_folder,'clone_positions_large_hand.dat']) == 2
    clone_positions{n_file} = load([main_folder_path{n_file},'/clone_positions.dat']);
    clone_positions_large{n_file} = load([an_folder,'clone_positions_large_hand.dat']);
    clone_positions_large_all_nodes{n_file} = load([an_folder,'clone_positions_large_all_nodes_hand.dat']);
elseif exist([main_folder_path{n_file},'/clone_positions.dat']) == 2
    clone_positions{n_file} = load([main_folder_path{n_file},'/clone_positions.dat']);
    clone_positions_large{n_file} = load([an_folder,'clone_positions_large.dat']);
    disp('No clone_positions_large_hand file found')
else
    clone_positions{n_file} = [];
    clone_positions_large{n_file} =[];
end

% duct width file
width_file_path = [an_folder,'/','edge_list_width_large.dat'];
node_width_file_path = [an_folder,'/','node_positions_width_large.dat'];
if isfile(width_file_path)
    edge_list_width_large{n_file} = load(width_file_path);
    node_positions_width_large{n_file} = load([an_folder,'node_positions_width_large.dat']);
else
    disp('Width datafile not found')
end

% lobe outline file
outline_file_path = [main_folder_path{n_file},'/','lobe_outline.dat'];
if isfile(outline_file_path)
    lobe_outline{n_file} = load(outline_file_path);
else
    disp('lobe outline datafile not found')
end

end

% REMOVE CLONES THAT CONTAIN ONLY DUCTAL CELLS AND MULTIPOTENT CLONES THAT
% SPAN A SINGLE DUCT
% clones will be removed if they only contain ductal cells OR if cells live
% in <= than min_number_ducts ducts.
min_number_ducts = 1;
 
if remove_ductal_clones == 1
    for n_file = 1:n_files
        if numel(clone_positions{n_file}) > 0
            [paths,dists,degrees] = node_level(edge_list_large{n_file},node_positions_large{n_file},source_node);
%             (clone_positions,min_number_ducts,dists)
            [clone_positions_temp,tot_num_clones,num_del_clones] = remove_ductal_only_clones_v1(clone_positions_large{n_file}, min_number_ducts,n_cell_types(n_file),dists);
%             [clone_positions_temp,tot_num_clones,num_del_clones] = remove_ductal_only_clones(clone_positions_large{n_file}, min_number_ducts,n_cell_types(n_file),dists,degrees,paths);
            [tot_num_clones,tot_num_clones-num_del_clones]
            cn = clone_positions_temp(:,1); cp = clone_positions_temp(:,7); ct = clone_positions_temp(:,6);
            
            cpos = clone_positions_temp(:,3:5);
%             cn(cp>max(node_positions_large{n_file}(:,1))) = []; cp(cp>max(node_positions_large{n_file}(:,1))) = [];
            clone_nodes_clean{n_file} = [cp, cn, ct, cpos];
            
            [clone_positions_duct,clone_positions_aci,clone_positions_myo] = split_ductal_acinar_cells(clone_positions_temp,n_cell_types(n_file));
            % duct
            cn = clone_positions_duct(:,1); cp = clone_positions_duct(:,7);
%             cn(cp>max(node_positions_large{n_file}(:,1))) = []; cp(cp>max(node_positions_large{n_file}(:,1))) = [];
            clone_nodes_clean_duct{n_file} = [cp, cn];
            % aci
            cn = clone_positions_aci(:,1); cp = clone_positions_aci(:,7);
%             cn(cp>max(node_positions_large{n_file}(:,1))) = []; cp(cp>max(node_positions_large{n_file}(:,1))) = [];
            clone_nodes_clean_aci{n_file} = [cp, cn];
            % myo
            if n_cell_types(n_file) > 3
            cn = clone_positions_myo(:,1); cp = clone_positions_myo(:,7);
%             cn(cp>max(node_positions_large{n_file}(:,1))) = []; cp(cp>max(node_positions_large{n_file}(:,1))) = [];
            clone_nodes_clean_myo{n_file} = [cp, cn];
            end
        else
            clone_nodes_clean{n_file} = [];
            clone_nodes_clean_duct = [];
            clone_nodes_clean_aci = [];
            clone_nodes_clean_myo = [];
        end
    end
end

%
% figure;
% for n_file = 1:n_files
%     subplot(n_files,1,n_file)
%     clone_properties(clone_positions{n_file},n_cell_types{n_file})
%     title(leg{n_file})
% end

%% plot reduced

for n_file = 1:n_files
    [~,dists,D] = node_level(edge_list_large{n_file},node_positions_large{n_file},source_node);
    
%     subplot(1,n_files,n_file)
    figure;
    h = plot_graph(edge_list_large{n_file},node_positions_large{n_file},source_node,0);
    title({leg{n_file},['#ducts: ',num2str(size(edge_list_large{n_file},1))]})
%     box off
    axis tight
    add_generations_to_tree(max(dists));
    set(gcf,'color','w');
    
    nbifs(n_file) = sum(D == 3);
end
% highlight(h,in,'Marker','s','NodeColor','b')

%% plot 3d
level_cutoff = 6;
plot_subtrees = 1;
dim = 3;
for n_file = 1:n_files
    figure;
[vtot(n_file),atot(n_file),vrudi(n_file),arudi(n_file),nodes_under_cutoff] = plot_3d_network_individual_branches(edge_list_large{n_file},node_positions_large{n_file},level_cutoff,source_node,plot_subtrees,dr(n_file,:),dim);
view(2) 
 
a = colorbar;
ylabel(a,'Increasing subtree size \rightarrow','FontSize',12,'Rotation',270);
set(a,'YTick',[]);
a.Label.Position(1) = 3;
a.Label.Position(1) = 2;
a.Label.Position(1) = 2.1;
legend('Rudimentary tree'); legend boxoff
xlabel('\mum'); ylabel('\mum'); zlabel('\mum')
% axis tight
set(gcf,'color','w');
axis equal

figure;
plot_tree_color_subtrees(edge_list_large{n_file},node_positions_large{n_file},level_cutoff,source_node,plot_subtrees)
a = colorbar;
ylabel(a,'Increasing subtree size \rightarrow','FontSize',12,'Rotation',270);
end

%% PLOT FULL 3D NETWORK WITH CLONES

for n_file = 1:n_files
    % asign random colours
    cn = clone_positions{n_file}(:,1);
    c = colormap('parula'); c = c(1:3:end,:);
    rp = randperm(size(c,1)); rp = rp(1:max(cn));
    c = c(rp,:);
   
    figure;
    % plot 3d network
    plot_3d_network(edge_list_large_all_nodes{n_file},node_positions_large_all_nodes{n_file},dr(n_file,:))
    hold on
    % add clones
    clone_nums = unique(clone_nodes_clean{n_file}(:,2));
    for i = 1:numel(clone_nums)
        in = find(clone_nodes_clean{n_file}(:,2) == clone_nums(i));
        c(i,:) = rand(1,3);
        scatter3(clone_nodes_clean{n_file}(in,4),clone_nodes_clean{n_file}(in,5),clone_nodes_clean{n_file}(in,6),20,'filled','MarkerFaceColor',c(i,:),'MarkerEdgeColor','k')
    end
    hold off
    axis equal
    box on
    xlabel('\mum');ylabel('\mum');ylabel('\mum')
    title(leg{n_file})
end

%% PLOT BRANCHING TREE WITH CLONES
clone_nodes = {};
for n_file = 1:n_files
    [~,dists,D] = node_level(edge_list_large{n_file},node_positions_large{n_file},source_node);
    
    cn = clone_nodes_clean{n_file}(:,2); cp = clone_nodes_clean{n_file}(:,1);
    clone_nodes{n_file} = [cp, cn];
    
    c = colormap('parula'); c = c(1:3:end,:);
    rp = randperm(size(c,1)); rp = rp(1:max(clone_nodes{n_file}(:,2)));
    c = c(rp,:);

    figure;
    % plot brnaching tree
    plot_graph_with_clones(edge_list_large{n_file},node_positions_large{n_file},source_node,0,clone_nodes{n_file},[],dr(n_file,:),c);
    % add level lines
    add_generations_to_tree(max(dists));
    title({leg{n_file},['#ducts: ',num2str(size(edge_list_large{n_file},1))]})
end

%% FUNCTIONS %%

function h = plot_graph(edge_list,node_positions,source_node,use_weights,weights)
m = max(node_positions(:,1));
edge_list(any(isnan(edge_list), 2), :) = m + 1;
% node_positions(any(isnan(node_positions), 2), :) = [];
if use_weights == 0
    G = graph(edge_list(:,1),edge_list(:,2));
    h = plot(G,'layout','layered','NodeLabel',{},'Marker','none','EdgeColor','k');
    layout(h,'layered','sources',source_node);
else
    G = graph(edge_list(:,1),edge_list(:,2),weights);
    weights = G.Edges.Weight;
    weights(~isfinite(weights)) = 0.1;
    weights(weights == 0) = 0.1;
    weights = log(weights); weights(weights < 0) = 10^-10;
    h = plot(G,'.','layout','layered','LineWidth',weights,'NodeLabel',{});
    layout(h,'layered','sources',source_node);
end
end % end plot_graph

function [paths,dists,D,G] = node_level(edge_list,node_positions,source_node) 
% calculate the level of each node and if such node is a branching or end point
% output:
% paths : paths of nodes from the shourse to every other node
% dists : distance (in node jumps) from the shourse to every other node
% D : degrees of each node
% G : graph object
G = graph(edge_list(:,1),edge_list(:,2));
paths = {};
dists = [];
for i = 1:size(node_positions,1)
    nod = node_positions(i,1);
    if nod~=i
        disp('Warning: incompatible node_positions file')
        break;
    end
    if isfinite(nod)
    [P,d] = shortestpath(G,source_node,nod);
    paths{i} = P;
    dists(i) = d;
    else
        paths{i} = nan;
        dists(i) = nan;
    end
end
D = degree(G);
end  % end node_level

function c = plot_graph_with_clones(edge_list,node_positions,source_node,d3,clone_nodes,clone_positions,dr,c)
m = max(node_positions(:,1));
edge_list(any(isnan(edge_list), 2), :) = m + 1;
if d3 == 1
    G = graph(edge_list(:,1),edge_list(:,2),edge_list(:,3),'OmitSelfLoops');
    
    x = node_positions(:,2).*dr(1);
    y = node_positions(:,3).*dr(2);
    z = node_positions(:,4).*dr(3);
    
    h = plot(G,'XData',x,'YData',y,'ZData',z);
    hold on
    clone_nums = unique(clone_positions(:,1));
    for i = 1:numel(clone_nums)
        in = find(clone_positions(:,1) == clone_nums(i));
        c(i,:) = rand(1,3);
        scatter3(clone_positions(in,3),clone_positions(in,4),clone_positions(in,5),50,'filled','MarkerFaceColor',c(i,:),'MarkerEdgeColor','k')
    end
    hold off
else
    G = graph(edge_list(:,1),edge_list(:,2),edge_list(:,3),'OmitSelfLoops');
    h = plot(G,'layout','layered','NodeLabel',{},'Marker','none','EdgeColor','k');
    layout(h,'layered','sources',source_node);
    un = unique(clone_nodes(:,2));
    hold on
    for i = 1:numel(un)
        ind = find(clone_nodes(:,2) == un(i));
        highlight(h,clone_nodes(ind,1),'Marker','o','NodeColor',c(i,:),'MarkerSize',3);
    end
    hold off
    axis tight
end

end % end plot_graph_with_clones

function add_generations_to_tree(max_level)
grid on
yticks(1:(max_level+1))
yTickLabels = arrayfun(@num2str,max_level:-1:0,'uni',false);
yticklabels(yTickLabels)

box off
ax1 = gca;  
ax1.XAxis.Visible = 'off';
% bg = get(gcf, 'Color');
% axes('xcolor', bg);
ylabel('Tree level (i)')
end % end add_generations_to_tree

function [clone_positions,tot_num_clones,num_del_clones] = remove_ductal_only_clones_v1(clone_positions,min_number_ducts,n_cell_types,dists)
% removes clones that are compodes exclusively of ductal cells
clone_index = clone_positions(:,1);
cell_type   = clone_positions(:,6);
clone_nodes = clone_positions(:,7);

indc_to_delete = [];
% run through all clones
un = unique(clone_index);
tot_num_clones = numel(un);
num_del_clones = 0;
for i = 1:numel(un)
    % find all cell i clone
    indc = find(clone_index == un(i));
    % extact cell types
    cell_type_in_clone = cell_type(indc);
    % find nodes/ducts of cell in clone
    nodes_in_clone = clone_nodes(indc); un_ducts = length(unique(nodes_in_clone));
    min_level = min(dists(nodes_in_clone));
    if n_cell_types == 3
        % acinar distal     1
        % acinar proximal	2
        % ductal            3
        if un_ducts <= min_number_ducts || min_level <= 1
            % store indices to delete
           indc_to_delete = [indc_to_delete; indc];
           num_del_clones = num_del_clones+1;
        end
    elseif n_cell_types == 8
        % acinar high Mist1	1
        % acinar low Mist1	2
        % ductal intercallated	3
        % ductal intralobular luminal	4
        % ductal intralobular basal	5
        % ductal interlobular luminal	6
        % ductal inerloboular basal	7
        % myoepithelial	8
        if un_ducts <= min_number_ducts || min_level <= 2
           indc_to_delete = [indc_to_delete; indc];
           num_del_clones = num_del_clones+1;
        end
    else
        disp(['Inconsistent number of cell types:',num2str(n_cell_types)]);
        break;
    end
    
end
% remove cells 
clone_positions(indc_to_delete,:) = [];
end % end remove_ductal_only_clones_v1


function [clone_positions_duct,clone_positions_aci,clone_positions_myo] = split_ductal_acinar_cells(clone_positions,n_cell_types)
clone_positions_duct = [];
clone_positions_aci = [];
clone_positions_myo = [];

cell_type   = clone_positions(:,6);
% n_cell_types = max(unique(cell_type)); % number of cell types

if n_cell_types == 3
    % acinar distal     1
    % acinar proximal	2
    % ductal            3
    ind_aci = find(cell_type <= 2);
    ind_duct = find(cell_type == 3);
    clone_positions_duct = clone_positions(ind_duct,:); 
    clone_positions_aci = clone_positions(ind_aci,:);
    clone_positions_myo = [];
elseif n_cell_types == 8
    % acinar high Mist1	1
    % acinar low Mist1	2
    % ductal intercallated	3
    % ductal intralobular luminal	4
    % ductal intralobular basal	5
    % ductal interlobular luminal	6
    % ductal inerloboular basal	7
    % myoepithelial	8
    ind_aci = find(cell_type <= 2);
    ind_duct = find(cell_type >= 3 & cell_type < 8);
    ind_myo = find(cell_type == 8);
    clone_positions_duct = clone_positions(ind_duct,:); 
    clone_positions_aci = clone_positions(ind_aci,:);
    clone_positions_myo = clone_positions(ind_myo,:);
else
end
end % end split_ductal_acinar_cells

function plot_3d_network(edge_list,node_positions,dr)
hold on
for i=1:length(edge_list(:,1))
    n1 = edge_list(i,1);
    n2 = edge_list(i,2);
    r1 = node_positions(n1,2:4).*dr;
    r2 = node_positions(n2,2:4).*dr;
%     plot([r1(1) r2(1)],[r1(2) r2(2)],'-k','linewidth',1.5)
    plot3([r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)],'-k')
end
hold off
end % end plot_3d_network

function [vtot,atot,vrudi,arudi,nodes_under_cutoff] = plot_3d_network_individual_branches(edge_list,node_positions,level_cutoff,source_node,plot_subtrees,dr,dim)
[~,dists,~,G] = node_level(edge_list,node_positions,source_node);

nodes_under_cutoff = find(dists <= level_cutoff);
nodes_in_cutoff = find(dists == level_cutoff);

cmap = colormap(parula((level_cutoff)^2));

for i = 2:4
    pos = node_positions(:,i);
    minp = min(pos);
    node_positions(:,i) = node_positions(:,i) - minp;
end

hold on
for i=1:length(edge_list(:,1))
    n1 = edge_list(i,1);
    n2 = edge_list(i,2);
    if all(ismember([n1,n2],nodes_under_cutoff))
        r1 = node_positions(n1,2:4);
        r2 = node_positions(n2,2:4);
        if dim == 2
            plot([r1(1) r2(1)].*dr(1),[r1(2) r2(2)].*dr(2),'-k','LineWidth',0.5)
        else
            plot3([r1(1) r2(1)].*dr(1),[r1(2) r2(2)].*dr(2),[r1(3) r2(3)],'-k','LineWidth',0.5)
        end
    end
end
rudnodes = nodes_under_cutoff; rudnodes(rudnodes == 1) = [];

[~,arudi] = boundary([node_positions(rudnodes,2).*dr(1),node_positions(rudnodes,3).*dr(2)]);%,node_positions(rudnodes,4)]);%
[~,atot] = boundary([node_positions(2:end,2).*dr(1),node_positions(2:end,3).*dr(2)]);%,node_positions(2:end,4)]); %
[~,vrudi] = boundary([node_positions(rudnodes,2).*dr(1),node_positions(rudnodes,3).*dr(2),node_positions(rudnodes,4)]);%
[~,vtot] = boundary([node_positions(2:end,2).*dr(1),node_positions(2:end,3).*dr(2),node_positions(2:end,4)]); %

if plot_subtrees == 1
    for nin = 1:length(nodes_in_cutoff)
        node = nodes_in_cutoff(nin);
        nodes_in_subtree = find_nodes_per_duct(dists,G,node,level_cutoff);
      
        nducts = 0;
        for i=1:length(edge_list(:,1))
            n1 = edge_list(i,1);
            n2 = edge_list(i,2);
            if all(ismember([n1,n2],nodes_in_subtree))
                r1 = node_positions(n1,2:4);
                r2 = node_positions(n2,2:4);
                %     plot([r1(1) r2(1)],[r1(2) r2(2)],'-k','linewidth',1.5)
                nducts = nducts + 1;
            end

        end
        subtree_size(nin) = nducts;
    end
    [~,I] = sort(subtree_size,'Descend');
   
    for n = 1:length(nodes_in_cutoff)
        nin = I(n);
        node = nodes_in_cutoff(nin);
        nodes_in_subtree = find_nodes_per_duct(dists,G,node,level_cutoff);
        
        nducts = 0;
        for i=1:length(edge_list(:,1))
            n1 = edge_list(i,1);
            n2 = edge_list(i,2);
            if all(ismember([n1,n2],nodes_in_subtree))
                r1 = node_positions(n1,2:4);
                r2 = node_positions(n2,2:4);
                if dim == 2
                    h(i) = plot([r1(1) r2(1)].*dr(1),[r1(2) r2(2)].*dr(2),'-','LineWidth',0.5,'Color',cmap(n,:));    
                else
                    h(i) = plot3([r1(1) r2(1)].*dr(1),[r1(2) r2(2)].*dr(2),[r1(3) r2(3)],'-','LineWidth',0.5,'Color',cmap(n,:));
                end
                nducts = nducts + 1;
            end
            subtree_size(i) = nducts;
        end
    end
end
hold off
end % end plot_3d_network_individual_branches

function plot_tree_color_subtrees(edge_list,node_positions,level_cutoff,source_node,plot_subtrees)
[~,dists,D,G] = node_level(edge_list,node_positions,source_node);

h = plot_graph(edge_list,node_positions,source_node,0);
axis tight
add_generations_to_tree(max(dists));
set(gcf,'color','w');

nodes_under_cutoff = find(dists <= level_cutoff);
nodes_in_cutoff = find(dists == level_cutoff);

cmap = colormap(parula((level_cutoff)^2));

if plot_subtrees == 1
    for nin = 1:length(nodes_in_cutoff)
        node = nodes_in_cutoff(nin);
        nodes_in_subtree = find_nodes_per_duct(dists,G,node,level_cutoff);
      
        nducts = 0;
        for i=1:length(edge_list(:,1))
            n1 = edge_list(i,1);
            n2 = edge_list(i,2);
            if all(ismember([n1,n2],nodes_in_subtree))
                r1 = node_positions(n1,2:4);
                r2 = node_positions(n2,2:4);
                nducts = nducts + 1;
            end

        end
        subtree_size(nin) = nducts;
    end
    [~,I] = sort(subtree_size,'Descend');
   
    for n = 1:length(nodes_in_cutoff)
        nin = I(n);
        node = nodes_in_cutoff(nin);
        nodes_in_subtree = find_nodes_per_duct(dists,G,node,level_cutoff);
        for i=1:length(edge_list(:,1))
            n1 = edge_list(i,1);
            n2 = edge_list(i,2);
            if all(ismember([n1,n2],nodes_in_subtree))
                highlight(h,[n1,n2],'EdgeColor',cmap(n,:));%,'LineWidth',1.5)  
            end
        end
    end
end
hold off
end % end plot_tree_color_subtrees
