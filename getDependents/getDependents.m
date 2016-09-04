classdef getDependents < handle
%GETDEPENDENTS Analyzes a file for dependent functions.
% This function will analyze an input file for all dependent functions.
% It works with scripts, functions, and classDef files.  It also can
% return any calls to external functions, toolboxes, built-in functions,
% Java functions, and other unknown functions.  The function was developed
% specifically to handle large object-oriented projects with many class
% definitions, packages, and associated functions.
%
% Included are methods to find all of the parent functions that call a
% listed function (or class).  It can optionally return  a sparse matrix
% indicating who calls whom.
%
% Finally, a method is provided to find all of the orphans of a folder or
% class package that are not called by anyone -- so called "orphans".
%
% For use with classes, this function must be called from the parent directory
% of any package folders, and other folders containing potential files of interest
% must be on the Matlab path.
%
% Example usage:
%    % Constructor
%    D = getDependents('Main.m');
%
%    % Return list of all dependent files
%    w = D.dependentFiles
%
%    % Return direct calls made by ProgramManager.m
%    x = D.whoAreCalledBy('+myPkg\ProgramManager.m')
%
%    % Return files that directly call Job.m
%    y = D.whoCall('+myPkg\Job.m')
%
%    % Return orphans in (or below) the +myPkg package
%    z = D.findOrp('+myPkg')
%
%    % Return orphans in (or below) the \utilities folder
%    z = D.findOrp('utilities')
%
%    % List external, toolbox, builtin, Java, and unknown calls
%    [ex,tb,bi,ja,un] = D.parseCalls('javaProgram.m')
%
% This function will recursively analyze any external functions called by
% the input rootFile.  It does not recursively search toolbox functions,
% built-in functions, Java functions, or other unknown functions.  It does
% not list other functions or methods within the same input rootFile.
%
% getDependents.m uses undocumented features of Matlab (mtree, mtfind)
% which may not be available in future releases.  Tested successfully with
% R2010a, R2011a, R2012a, and R2013a.
%
% Requires "Recursive Directory Listing - Enhanced RDIR" by Thomas Vanaret
% available on Matlab FileExchange.

% Version: 1.0
% Author:  Mark W. Brown
% Date:    04 Sept 2013

    properties
        dependentFiles = {};   % list of all M-files required by rootFile
        dependencyMatrix;      % dependencyMatrix(i,j)=1 if dependentFiles{i} calls dependentFiles{j}
    end
    
    methods
        function self = getDependents(rootFile)
            %GETDEPENDENTS Main constructor.
            % getDependents returns an object containing a list of all
            % dependent files (starting with the input rootFile), and a
            % matrix which shows how the files depend on each other.
            
            fprintf('Searching for dependent files...')
            % Initialize processed file list:
            processed = {};

            % Initialize the list of files to process:
            pending = {which(rootFile)};

            % Keep processing until there is nothing in the pending list:
            while ~isempty(pending)
                % Pick off the top of the pending list:
                currFile = pending{1};
                % If it hasn't already been processed, then parse it:
                if ~ismember(currFile,processed)
                    % Add current file to the processed list:
                    processed = [processed; currFile];
                    % Remove it from the pending list:
                    pending(1) = [];
                    % Find the new files that currFile depends on:
                    newFiles = self.parseCalls(currFile);
                    % And add them to the pending list:
                    pending = [pending; newFiles];
                else
                    % If it has been processed, remove it from pending:
                    pending(1) = [];
                end
            end
            
            % Sort (in Windows Explorer sort order) and return the results:
            [~,ind] = sort(lower(processed));
            self.dependentFiles = processed(ind);
            fprintf('  Finished.\n')
            
            % Call traceAll to fill out the dependencyMatrix
            self.traceAll;
        end
        
        function traceAll(self)
            %TRACEALL This method builds the matrix of direct calls.
            % The entry (i,j)denotes that file{i} directly calls file{j}.
            fprintf('Tracing direct calls...')
            
            % Get the list of all dependent files:
            infiles = self.dependentFiles;
            numFiles = length(infiles);
            
            % Initialize a sparse matrix of file dependencies:
            S = sparse(numFiles, numFiles);

            % Loop through all the dependent files:
            for ii = 1:numFiles
                % For each file:
                currFile = infiles{ii};
                % Every file depends on itself:
                % S(ii,ii) = 1;
                % Find the files it calls directly:
                newFiles = self.parseCalls(currFile);
                % Loop through those references:
                for jj = 1:length(newFiles)
                    % Find it's position in the matrix:
                    ind = strcmp(newFiles{jj},infiles);
                    % And record it:
                    S(ii,ind) = 1;
                end
            end
            
            % Record the dependency matrix:
            self.dependencyMatrix = S;
            fprintf('  Finished.\n')
        end
        
        function C = whoAreCalledBy(self,filespec)
            %WHOARECALLEDBY Returns a list of files called by the input file.
            
            % Expand to a full filespec:
            filespec = which(filespec);
            
            % Retrieve the file dependency matrix:
            S = self.dependencyMatrix;

            % Find the row position of filespec in the file list:
            infiles = self.dependentFiles;
            ind = strcmp(filespec,infiles);

            % Get the indices of files called by filespec:
            calls = find(S(ind,:));
            
            if nargout
                % Return the results in a cell array:
                C = cell(length(calls),1);
                for ii = 1:length(calls)
                    C{ii} = infiles{calls(ii)};
                end
            else
                % Display the results to the command window:
                for ii = 1:length(calls)
                    fprintf('    %s\n',infiles{calls(ii)});
                end
            end
        end
        
        function C = whoCall(self,filespec)
            %WHOCALL Returns a list of files that call the input file.
            
            % Expand to a full filespec:
            filespec = which(filespec);
            
            % Retrieve the file dependency matrix:
            S = self.dependencyMatrix;

            % Find the column position of filespec in the file list:
            infiles = self.dependentFiles;
            ind = strcmp(filespec,infiles);

            % Get the indices of files called by filespec:
            calls = find(S(:,ind));
            
            if nargout
                % Return the results in a cell array:
                C = cell(length(calls),1);
                for ii = 1:length(calls)
                    C{ii} = infiles{calls(ii)};
                end
            else
                % Display the results to the command window:
                for ii = 1:length(calls)
                    fprintf('    %s\n',infiles{calls(ii)});
                end
            end
        end
        
        function orp = findOrp(self, rootFolder)
            %FINDORP Finds all of the apparent orphans in (or below) the root folder.
            % Orphans are files that are not directly
            % called by any of the dependent files.
            %
            % Inputs:
            %    rootFolder == folder to begin search for orphans
            %    (defaults to the current folder if not specified).
            %
            % Outputs:
            %    orp == list of apparent orphans (optional)
            %
            % Requires:  rdir.m (recursive directory listing)
            
            % Determine where to start:
            if nargin>1
                temp = what(rootFolder);
                startFolder = temp.path;
            else
                startFolder = cd;
            end

            % Get list of all files at or below the current folder:
            dirStruct = rdir([startFolder,'\**\*.m']);
            fileCell = struct2cell(dirStruct);
            allFiles = fileCell(1,:)';

            if nargout
                % Initialize the list of orphans:
                orp = {};

                % Loop thru all of the files and check to see if they
                % are on the dependent files list.  If so, skip them.
                % Otherwise, put them on the orphan list.
                for ii = 1:length(allFiles)
                    file = allFiles{ii};
                    if ~ismember(file,self.dependentFiles);
                        orp{length(orp)+1} = file;
                    end
                end

                % Sort and transpose the results:
                [~,ind] = sort(lower(orp));
                orp = orp(ind)';
            else
                % Loop thru all of the files and check to see if they
                % are on the dependent files list.  If so, skip them.
                % Otherwise, display them.
                for ii = 1:length(allFiles)
                    file = allFiles{ii};
                    if ~ismember(file,self.dependentFiles);
                        fprintf('%s\n',file);
                    end
                end
            end
        end


        function [extFiles, tbFiles, biFuncs, jaFuncs, unkFuncs] = parseCalls(self,filename)
        %PARSECALLS Parses the input file for function calls.
        % This function is adapted from Matlab's getcallinfo.m.  If
        % filename calls overloaded methods, this function may not return
        % the right results.  The method returned is the same as would be
        % obtained by typing >> which(<methodname>) at the command line.
        % The function returns a list of all external files (extFiles),
        % Matlab Toolbox files (tbFiles), Matlab builtin functions
        % (biFuncs), Java functions (jaFuncs), and unknown functions
        % (unkFuncs).

        % Parse m-file into tree structure:
        tree = mtree(filename,'-file');

        % Find calls
        calls = mtfind(tree,'Kind','CALL');

        % Get normal function calls:
        callInfo = calls.Left.strings;

        % Get the potential dot calls:
        dots = mtfind(tree,'Kind','ID', 'Parent.Kind', 'DOT');
        idx = dots.indices;
        dotCalls = cell(1,length(idx));
        for i=1:length(idx)
            top = tree.select(idx(i));
            dotName = string( top );
            while (iskind( Parent(top), 'DOT'))
                top = Parent(top);
                dotName = [dotName '.' string( Right(top))];
            end
            dotCalls{i} = dotName; 
        end
        callInfo = [callInfo, dotCalls];

        % Get superclass references:
        superClass = tree.Cexpr.Right;
        callInfo = [callInfo, superClass.strings];

        % Get calls to other constructors:
        otherConst = calls.Left.Right;
        callInfo = [callInfo, otherConst.strings];

        % Get embedded function references (@callbacks):
        ATs = mtfind(tree,'Kind','AT');
        funcRef = ATs.Arg;
        callInfo = [callInfo, funcRef.strings];

        % Get unique results (and transpose):
        callInfo = unique(callInfo');

        % Initialize output variables:
        extFiles = {};  % external files
        tbFiles = {};   % toolbox files
        biFuncs = {};   % built-in functions
        jaFuncs = {};   % java functions
        unkFuncs = {};  % unknown functions

        % Get ML Root:
        r = matlabroot;

        % Get names of files/functions called:
        for ii = 1:length(callInfo)
            w = which(callInfo{ii});
            if ~isempty(w)
                if strcmp(w(1:4),'java')
                    % Java functions
                    jaFuncs{length(jaFuncs)+1} = w;
                elseif strcmp(w(1:8),'built-in')
                    % Built-in functions
                    if strfind(w,'undocumented')
                         biFuncs{length(biFuncs)+1} = [w,' ',callInfo{ii}];
                    else
                        biFuncs{length(biFuncs)+1} = w;
                    end
                elseif strcmp(w(1:min(length(r),length(w))),r(1:min(length(r),length(w))))
                    % Toolbox m-files
                    tbFiles{length(tbFiles)+1} = w;
                elseif exist(w) == 2
                    % External (user) m-files
                    if ~strcmp(w,filename) % don't include self
                        extFiles{length(extFiles)+1} = w;
                    end
                else
                    % Unknown functions:
                    unkFuncs{length(unkFuncs)+1} = w;
                end
            end
        end

        % Transpose and return:
        extFiles = unique(extFiles');
        tbFiles = unique(tbFiles');
        biFuncs = unique(biFuncs');
        jaFuncs = unique(jaFuncs');
        unkFuncs = unique(unkFuncs');

        end
    end
end
