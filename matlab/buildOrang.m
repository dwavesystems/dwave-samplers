function buildOrang( boostDir, extraFlags )

clear orang_minsum orang_sample orang_greedyvarorder

if nargin >= 1 && ~isempty(boostDir)
  if ~ischar(boostDir) || size(boostDir, 1) ~= 1 || ~exist(boostDir, 'dir')
    error('boostDir must be an existing directory')
  end
end

if nargin >= 2
  if ischar(extraFlags)
    extraFlags = { extraFlags };
  end

  if ~iscellstr(extraFlags)
    error('extraFlags must be a string or cell array of strings')
  end
end

srcFileSets = {
  {'orang_minsum.cpp' 'orang_mex.cpp'} ...
  {'orang_mincount.cpp' 'orang_mex.cpp'} ...
  {'orang_sample.cpp' 'orang_mex.cpp'} ...
  {'orang_greedyvarorder.cpp' 'orang_mex.cpp'} ...
};

outDir = fileparts(mfilename('fullpath'));
orangIncl = [ '-I' fullfile(fileparts(outDir), 'src') ];

mexArgs = { '-O', '-largeArrayDims', '-outdir', outDir, orangIncl };
if nargin >= 1 && ~isempty(boostDir)
  mexArgs = [ mexArgs [ '-I' boostDir ] ];
end
if nargin >= 2
  mexArgs = [ mexArgs extraFlags(:)' ];
end

for ii=1:numel(srcFileSets)
  srcFiles = cellfun(@(s) fullfile(outDir, s), srcFileSets{ii}, 'UniformOutput', false);
  [zzz basename] = fileparts(srcFiles{1});
  fprintf('Building %s...\n', basename);

  mex(mexArgs{:}, srcFiles{:})
end


