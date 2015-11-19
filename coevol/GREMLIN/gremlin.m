function model = gremlin(fname_train, outfile, varargin);
if(nargin>2)
  opts = varargin;
else
  opts = cell(0);
end;

if(mod(length(opts),2)~=0)
  error('Incorrect number of input parameters!');
  return;
end;
options = parse_gremlin_options(opts)

X  = read_msa(fname_train);
nStates= max(max(X));
if(isfield(options, 'reweight'))
	options.seqDepWeights = (1./(1+sum(squareform(pdist(X, 'hamm')<options.reweight))))';
  disp(['reweighing sequences n ' num2str(size(X,1)) ' neff ' num2str(sum(options.seqDepWeights)) ]);
end
if (isfield(options, 'weight_file'))
  options.seqDepWeights = load(options.weight_file);
end;

start_time=tic;
if(isfield(options, 'lambda_mat'))
  lambda_vector = squareform(options.lambda_mat-diag(diag(options.lambda_mat)));
  options.lambda_vector=lambda_vector;
  disp(['initialized lambda_vector from lambda mat ' num2str(size(options.lambda_vector))]);
end;

model = LLM2_train(X, options);

toc(start_time)

[~, ~, edge_norm2_mat, ~] = compute_edge_norms(model, max(max(X)), size(X,2),outfile);

if(options.apc==1)
  edge_norm2_mat = apc(edge_norm2_mat);
end;

write_to_file(edge_norm2_mat, outfile,1)
