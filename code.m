	%group = 'L';
	
	groupvector = cell(4,1);
	groupvector{1} = 'H04';
	groupvector{2} = 'H10';
	groupvector{3} = 'M';
	groupvector{4} = 'L';
	
	for groupcounter = 1:4
	group = groupvector{groupcounter};
	
	%Generate file names
	Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
	
	filenames_letters = num2str(zeros(96,1));
	filenames_numbers = zeros(96,1);
	letter = 1;
	for i = 1:12:96
	filenames_letters(i:i+11) = Alphabet(letter);
	letter = letter+1;
	for j = 0:11
	filenames_numbers(i+j) = j+1;
	end
	end
	filenames_numbers = num2str(filenames_numbers);
	
	for i = 1:96
	if str2num(filenames_numbers(i,:)) < 10
	filenames_numbers(i,:) = ['0' filenames_numbers(i,2)];
	end
	end
	
	filenames = cellstr(num2str(zeros(96,1)));
	wellnames = filenames;
	c = ['_', group ,'_SAPA-23.csv'];
	cID = ['_', group ,'_SAPA-23'];
	
	for i = 1:length(filenames)
	plateID(i) = {[filenames_letters(i) filenames_numbers(i,:)]}';
	filenames(i) = {[[filenames_letters(i) filenames_numbers(i,:)], c]};
	wellnames(i) = {[[filenames_letters(i) filenames_numbers(i,:)], cID]};
	end
	
	% Rip sequences
	seq_long = cell(96,1);
	foldername = [group, '_sapa23_csv\'];
	for i = 1:96
	[num, txt] = xlsread([foldername, filenames{i}], 'A2:A2');
	seq_long(i) = txt;
	end
	
	% Define sequences
	seq_plasmid = [[See other appendices]];
	aa_plasmid = DNAtranslator(seq_plasmid(1:end-2));
	aa_wt = [[See other appendices]];
	seq_wt = seq_plasmid(668*3-2 : 668*3-2+117*3-1);
	aa_comp = [[See other appendices]];
	seq_comp = seq_plasmid(689*3-2 : 689*3-2+length(aa_comp)*3-1);
	
	preZ = seq_plasmid(2064-19:2064);
	postZ = seq_plasmid(2065+length(seq_comp):2065+length(seq_comp)+19);
	
	% Trim and classify
	seq = seq_long;
	comment = cell(96,1);
	error_fail = 'Sequencing failed';
	error_short = 'Sequence too short to align successfully';
	error_N = 'Z-domain contains undetermined bases. Check reverse sequence';
	error_indel = 'Sequence contains insertions or deletions';
	error_wt = 'Sequence is wild-type';
	error_align = 'Alignment failed';
	
	for i = 1:96
	% Cut-out comp. seq.
	if length(seq_long{i}) >= 174+2*20
	Zstart = strfind(seq_long{i}, preZ) + length(preZ);
	Zend = strfind(seq_long{i}, postZ) -1;
	seq{i} = seq_long{i}(Zstart:Zend);
	elseif length(seq_long{i}) == 5 && mean(seq_long{i} == 'NNNNN') == 1
	comment{i} = error_fail;
	seq{i} = seq_long{i};
	elseif isempty(comment{i}) == 1 && length(seq_long{i}) < 174+2*20
	comment{i} = error_short;
	seq{i} = seq_long{i};
	else
	comment{i} = error_align;
	seq{i} = seq_long{i};
	end
	% Is comp. seq. longer/shorter/wt?
	if isempty(comment{i}) == 1 && length(seq{i}) ~= length(seq_comp)
	comment{i} = error_indel;
	elseif isempty(comment{i}) == 1 && length(seq{i}) == length(seq_comp)
	if seq{i} == seq_comp
	comment{i} = error_wt;
	end
	end
	% Does comp. seq. contain bases other than ACTG?
	if isempty(comment{i}) == 1 && length(seq_long{i}) >= 174 ...
	&& length(strfind(seq{i}, 'A')) + length(strfind(seq{i}, 'T')) ...
	+ length(strfind(seq{i}, 'C')) + length(strfind(seq{i}, 'G')) < length(seq{i})
	comment{i} = error_N;
	elseif isempty(comment{i}) ~= 1 && length(seq_long{i}) >= 174 && ...
	length(strfind(seq{i}, 'A')) + length(strfind(seq{i}, 'T')) ...
	+ length(strfind(seq{i}, 'C')) + length(strfind(seq{i}, 'G')) < length(seq{i})
	comment{i} = [comment{i},' & ', error_N];
	end
	end
	
	% Translate sequences
	aa = cell(96,1);
	counter = 0;
	for i = 1:96
	counter = counter +1;
	if isempty(comment{i}) == 1
	aa{i} = DNAtranslator(seq{i});
	else
	aa{i} = '-';
	end
	end
	
	% Compare to wt seq
	subsIndex = cell(96,1);
	subsFrom = subsIndex;
	subsTo = subsIndex;
	
	for i = 1:96
	if isempty(comment{i}) == 1
	subsIndex{i} = find((aa{i} == aa_comp) == 0);
	subsFrom(i,1) = {aa_comp(subsIndex{i,1})};
	subsTo(i,1) = {aa{i}(subsIndex{i,1})};
	else
	subsFrom(i,1) = {'-'};
	subsTo(i,1) = {'-'};
	subsIndex{i,1} = '-';
	end
	end
	
	% Switch index to position
	position = cell(96,1);
	for i = 1:96
	if isnumeric(subsIndex{i}) == 1
	position{i} = subsIndex{i} + 21;
	else
	position{i} = subsIndex{i};
	end
	
	if length(position{i}) > 1
	position_merge = zeros(1,length(position{i}));
	for j = 1:length(position{i})
	position_merge(j) = position{i}(j);
	end
	position{i} = num2str(position_merge);
	comment{i} = 'Multiple substitutions';
	end
	end
	
	% Sort
	
	for i = 1:96
	if isnumeric(position{i}) == 1
	position_sort(i,1) = position{i};
	else
	position_sort(i,1) = 0;
	end
	end
	
	[sortVal, sortIndex] = sort(position_sort);
	
	output = cell(96,8);
	output(:,1) = wellnames(:,1);
	output(:,2) = comment(:,1);
	output(:,3) = seq_long(:,1);
	output(:,4) = seq(:,1);
	output(:,5) = aa(:,1);
	output(:,6) = position(:,1);
	output(:,7) = subsFrom(:,1);
	output(:,8) = subsTo(:,1);
	
	output_sorted = cell(96,8);
	
	for i = 1:96
	for j = 1:8
	output_sorted{i,j} = output{sortIndex(i),j};
	end    
	end
	
	headers = {'Well name', 'Comment', 'Full NT sequence' ...
	'Aligned NT sequence', 'Aligned AA sequence', 'AA position'...
	'From AA', 'To AA'};
	output_sorted
	
	xlswrite('results.xls', headers, group, 'A1:H1')
	xlswrite('results.xls', output_sorted, group, 'A2:H97')
	
	end