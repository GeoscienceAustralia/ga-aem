classdef gatdaem1d_system_handle < handle
    %GATDAEM1D_SYSTEM_HANDLE  MATLAB wrapper (handle class) for the gatdaem1d C API (gatdaem1d.dll/.so).
    %
    % This class wraps an opaque C/C++ system handle (uint64) returned by createhandle(),
    % and exposes convenience methods for common metadata queries and forward modelling.
    %
    % The wrapper uses MATLAB's loadlibrary/calllib mechanism, with prototypes/struct
    % definitions generated in gatdaem1d_proto.m from the public C header gatdaem1d.h.
    %
    % Basic usage
    % -----------
    %   S = gatdaem1d_system_handle('path/to/system.stm');
    %   nw = S.nwindows();
    %   [tLow,tHigh] = S.windowtimes();
    %
    %   % Geometry (must contain all fields of struct gatdaem1d_system_geometry)
    %   G = struct('tx_height',30,'tx_roll',0,'tx_pitch',0,'tx_yaw',0, ...
    %              'txrx_dx',0,'txrx_dy',0,'txrx_dz',0, ...
    %              'rx_roll',0,'rx_pitch',0,'rx_yaw',0);
    %
    %   % Earth (must contain E.conductivity [nlayersx1] and E.thickness [(nlayers-1)x1])
    %   E = struct('conductivity',[0.01; 0.001; 0.05], 'thickness',[20; 40]);
    %
    %   R = S.forwardmodel(G,E);    % returns struct with fields PX,PY,PZ,SX,SY,SZ
    %
    % Lifetime / cleanup
    % ------------------
    % This class does NOT implement a MATLAB delete() method. That means MATLAB will not
    % automatically call deletehandle() when the object is cleared. Call free_system_handle()
    % yourself when you are finished with the underlying C handle.
    %
    %   S = gatdaem1d_system_handle(systemfile);
    %   c = onCleanup(@() S.free_system_handle());
    %   ... use S ...
    %
    % ABI/version guarding
    % -------------------
    % On first load, ensureLibraryIsLoaded_() compares the ABI version reported by the
    % loaded library (gatdaem1d_get_abi_version) with the ABI constant embedded in the
    % generated proto (enuminfo.GATDAEM1D_EXPECTED_ABI...). If they differ, an error is thrown.
    %
    % Inputs/outputs (high level)
    % ---------------------------
    % Geometry G fields (all scalar, real):
    %   tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw
    %
    % Earth E fields:
    %   Required: conductivity (nlayers x 1), thickness ((nlayers-1) x 1)
    %   Optional IP: iptype plus (chargeability,timeconstant,frequencydependence), each (nlayers x 1)
    %
    % Notes
    % -----
    % * All numeric inputs are passed as double (geometry/earth arrays) or int32 (enums/counts).
    % * If a library call fails, throwIfNotOk_() raises an error with the status code and
    %   the library's last error string (gatdaem1d_get_last_error_string).
    %
    % MATLAB RAII wrapper for gatdaem1d.dll based on gatdaem1d_proto.m
    %
    % Functions wrappedfstatus

    %   gatdaem1d_get_last_error_string(int32)
    %   createhandle(cstring, uint64Ptr)
    %   deletehandle(uint64Ptr)
    %   nsamplesperwaveform(uint64, int32Ptr)
    %   waveform(uint64, doublePtr, doublePtr, doublePtr)
    %   nwindows(uint64, int32Ptr)
    %   nturns(uint64, int32Ptr)
    %   looparea(uint64, doublePtr)
    %   basefrequency(uint64, doublePtr)
    %   peakcurrent(uint64, doublePtr)
    %   windowtimes(uint64, doublePtr, doublePtr)
    %   forwardmodel(uint64, gatdaem1d_system_geometryPtr, gatdaem1d_earthPtr, gatdaem1d_responsePtr)
    %   derivative(uint64, gatdaem1d_system_geometryPtr, gatdaem1d_earthPtr, cstring, int32, gatdaem1d_responsePtr)
    %   fm_dlogc(uint64, gatdaem1d_system_geometryPtr, gatdaem1d_earthPtr, doublePtr)
    %   iptype_from_string(cstring, int32Ptr)

    properties (SetAccess = private)
        LibAlias (1,:) char = 'gatdaem1d'
        DllPath  (1,:) char = 'gatdaem1d.dll'
        ProtoFcn = @gatdaem1d_proto
        HandlePtr               % libpointer('uint64Ptr', uint64(0))
        SystemFile = [];
        Enums      = [];
        LibraryABIVersion = [];
        ProtoABIVersion   = [];
    end

    properties (Dependent)
        HandleValue uint64
        IsValid (1,1) logical
    end

    properties (Constant)

    end

    methods
        function abi_version = gatdaem1d_get_abi_version(obj)            
            obj.ensureObjectIsLoadedAndValid_();
            abi_version = calllib(obj.LibAlias, 'gatdaem1d_get_abi_version');            
        end
        function msg = gatdaem1d_get_last_error_string(obj, status)
            %GET_LAST_ERROR_STRING Return the library's last error message for a status code.
            %
            % C prototype:
            %   extern const char* __cdecl gatdaem1d_get_last_error_string(const gatdaem1d_status s);
            %
            % Thunk info (from your prototype file):
            %   thunkname: cstringint32Thunk
            %   LHS:       cstring
            %   RHS:       int32

            obj.ensureObjectIsLoadedAndValid_();
    
            if nargin < 2
                status = int32(obj.Enums.STATUS.GATDAEM1D_STATUS_OK);
            else
                status = int32(status);
            end

            % calllib returns a MATLAB char row vector for LHS='cstring'
            msg = calllib(obj.LibAlias, 'gatdaem1d_get_last_error_string', status);

            % Defensive cleanup
            if isempty(msg)
                msg = '';
            elseif isstring(msg)
                msg = char(msg);
            end
        end
        function obj = gatdaem1d_system_handle(systemfile, varargin)
            %GATDAEM1D_SYSTEM_HANDLE Construct and allocate a new underlying C handle.
            %
            %   S = gatdaem1d_system_handle(systemfile)
            %   S = gatdaem1d_system_handle(systemfile,'LibAlias',alias,'DllPath',dll,'ProtoFcn',@gatdaem1d_proto)
            %
            % Inputs
            %   systemfile : path to a system description file understood by the C++ library.
            %
            % Name-Value pairs
            %   LibAlias : library alias used by loadlibrary/calllib (default 'gatdaem1d').
            %   DllPath  : path to the shared library (default 'gatdaem1d.dll').
            %   ProtoFcn : function handle for the generated prototype (default @gatdaem1d_proto).
            %
            % Side effects
            %   - Loads the library if required (ensureLibraryIsLoaded_).
            %   - Calls createhandle(systemfile,&handle).
            %   - Throws on non-OK status or if a zero handle is returned.
            p = inputParser;
            p.addRequired('systemfile', @(s)ischar(s)||isstring(s));
            p.addParameter('LibAlias','gatdaem1d', @(s)ischar(s)||isstring(s));
            p.addParameter('DllPath','gatdaem1d.dll', @(s)ischar(s)||isstring(s));
            p.addParameter('ProtoFcn',@gatdaem1d_proto, @(f)isempty(f) || isa(f,'function_handle'));
            p.parse(systemfile, varargin{:});

            obj.LibAlias = char(p.Results.LibAlias);
            obj.DllPath  = char(p.Results.DllPath);
            obj.ProtoFcn = p.Results.ProtoFcn;

            obj.HandlePtr = libpointer('uint64Ptr', uint64(0));
            obj.SystemFile = systemfile;

            obj.ensureLibraryIsLoaded_();

            % Explicit cstring (avoid implicit marshaling)
            stmPtr = libpointer('cstring', char(systemfile));
            status = calllib(obj.LibAlias, 'createhandle', stmPtr, obj.HandlePtr);
            obj.throwIfNotOk_(status, 'createhandle');

            if obj.HandlePtr.Value == 0
                error('gatdaem1d_system_handle:CreateFailed', 'createhandle returned OK but handle is 0');
            end
        end
        function free_system_handle(obj)
            %FREE_SYSTEM_HANDLE Release the underlying C handle (deletehandle).
            %
            %   S.free_system_handle()
            %
            % This sets the C-side handle to 0 via deletehandle(&handle). Once freed,
            % IsValid will return false and most other methods will error.
            %
            % Recommended pattern:
            %   S = gatdaem1d_system_handle(systemfile);
            %   c = onCleanup(@() S.free_system_handle());
            try
                if ~isempty(obj.HandlePtr) && isa(obj.HandlePtr,'lib.pointer') && libisloaded(obj.LibAlias) && obj.HandlePtr.Value ~= 0
                    status = calllib(obj.LibAlias, 'deletehandle', obj.HandlePtr);                    
                    obj.throwIfNotOk_(status, 'nturns');
                end
            catch

            end
        end
        function v = get.HandleValue(obj)
            if isempty(obj.HandlePtr), v = uint64(0);
            else, v = uint64(obj.HandlePtr.Value);
            end
        end
        function tf = get.IsValid(obj)
            tf = ~isempty(obj.HandlePtr) && isa(obj.HandlePtr,'lib.pointer') && (obj.HandlePtr.Value ~= 0);
        end

        %% ---- error helpers ----
        function clear_last_error(obj)
            obj.ensureObjectIsLoadedAndValid_();
            calllib(obj.LibAlias, 'gatdaem1d_clear_last_exception_error');
        end
        function s = last_exception_string(obj)
            obj.ensureObjectIsLoadedAndValid_();
            s = calllib(obj.LibAlias, 'gatdaem1d_get_last_exception_error_string');
        end
        function s = status_string(obj, status)
            obj.ensureObjectIsLoadedAndValid_();
            s = calllib(obj.LibAlias, 'gatdaem1d_status_string', int32(status));
        end
        
        function iptype = iptype_from_string(obj, iptype_string)
            obj.ensureObjectIsLoadedAndValid_();
            if ~(ischar(iptype_string) || isstring(iptype_string))
                error('ip_type_from_string:BadType', 'iptype_string must be char or string.');
            end
            iptype_stringPtr = libpointer('cstring', char(iptype_string));
            iptypePtr        = libpointer('int32Ptr', int32(0));
            status = calllib(obj.LibAlias, 'iptype_from_string', iptype_stringPtr, iptypePtr);
            obj.throwIfNotOk_(status, 'iptype_from_string');
            iptype = int32(iptypePtr.Value);
        end

        %% ---- scalar queries ----
        function n = nsamplesperwaveform(obj)
            %NSAMPLESPERWAVEFORM Number of samples in the transmitter waveform arrays.
            %
            %   n = S.nsamplesperwaveform()
            %
            % This is used to size the buffers returned by waveform().
            obj.ensureObjectIsLoadedAndValid_();
            nPtr = libpointer('int32Ptr', int32(0));
            status = calllib(obj.LibAlias, 'nsamplesperwaveform', obj.HandleValue, nPtr);
            obj.throwIfNotOk_(status, 'nsamplesperwaveform');
            n = int32(nPtr.Value);
        end
        function n = nwindows(obj)
            %NWINDOWS Number of time windows in the system response.
            %
            %   nw = S.nwindows()
            obj.ensureObjectIsLoadedAndValid_();
            obj.ensureObjectIsLoadedAndValid_();
            nPtr = libpointer('int32Ptr', int32(0));
            status = calllib(obj.LibAlias, 'nwindows', obj.HandleValue, nPtr);
            obj.throwIfNotOk_(status, 'nwindows');
            n = int32(nPtr.Value);
        end
        function n = nturns(obj)
            %NTURNS Number of transmitter turns reported by the system definition.
            %
            %   n = S.nturns()
            obj.ensureObjectIsLoadedAndValid_();
            nPtr = libpointer('int32Ptr', int32(0));
            status = calllib(obj.LibAlias, 'nturns', obj.HandleValue, nPtr);
            obj.throwIfNotOk_(status, 'nturns');
            n = int32(nPtr.Value);
        end

        function a = looparea(obj)
            %LOOPAREA Transmitter loop area reported by the system definition.
            obj.ensureObjectIsLoadedAndValid_();
            aPtr = libpointer('doublePtr', 0);
            status = calllib(obj.LibAlias, 'looparea', obj.HandleValue, aPtr);
            obj.throwIfNotOk_(status, 'looparea');
            a = double(aPtr.Value);
        end

        function f = basefrequency(obj)
            %BASEFREQUENCY Base frequency reported by the system definition.
            obj.ensureObjectIsLoadedAndValid_();
            fPtr = libpointer('doublePtr', 0);
            status = calllib(obj.LibAlias, 'basefrequency', obj.HandleValue, fPtr);
            obj.throwIfNotOk_(status, 'basefrequency');
            f = double(fPtr.Value);
        end

        function i = peakcurrent(obj)
            %PEAKCURRENT Peak current reported by the system definition.
            obj.ensureObjectIsLoadedAndValid_();
            iPtr = libpointer('doublePtr', 0);
            status = calllib(obj.LibAlias, 'peakcurrent', obj.HandleValue, iPtr);
            obj.throwIfNotOk_(status, 'peakcurrent');
            i = double(iPtr.Value);
        end
        
        function [t, I, V] = waveform(obj)
            %WAVEFORM Get waveform sample vectors from the loaded system.
            %
            %   [t, I, V] = S.waveform()
            %
            % Outputs are column vectors of length nsamplesperwaveform():
            %   t : sample times
            %   I : current waveform
            %   V : voltage waveform
            obj.ensureObjectIsLoadedAndValid_();
            n = double(obj.nsamplesperwaveform());
            tPtr = libpointer('doublePtr', zeros(n,1));
            iPtr = libpointer('doublePtr', zeros(n,1));
            vPtr = libpointer('doublePtr', zeros(n,1));
            status = calllib(obj.LibAlias, 'waveform', obj.HandleValue, tPtr, iPtr, vPtr);
            obj.throwIfNotOk_(status, 'waveform');
            t = tPtr.Value; I = iPtr.Value; V = vPtr.Value;
        end

        function [tLow, tHigh] = windowtimes(obj)
            %WINDOWTIMES Get the low/high time bounds for each window.
            %
            %   [tLow,tHigh] = S.windowtimes()
            %
            % Outputs are column vectors of length nwindows().
            obj.ensureObjectIsLoadedAndValid_();
            n = double(obj.nwindows());
            lowPtr  = libpointer('doublePtr', zeros(n,1));
            highPtr = libpointer('doublePtr', zeros(n,1));
            status = calllib(obj.LibAlias, 'windowtimes', obj.HandleValue, lowPtr, highPtr);
            obj.throwIfNotOk_(status, 'windowtimes');
            tLow = lowPtr.Value; tHigh = highPtr.Value;
        end
        
        function R = forwardmodel(obj, G, E)
            %FORWARDMODEL Compute the forward response for geometry G and earth model E.
            %
            %   R = S.forwardmodel(G,E)
            %
            % Inputs
            %   G : struct with fields matching gatdaem1d_system_geometry (see class help).
            %   E : struct with required fields conductivity/thickness (and optional IP fields).
            %
            % Output
            %   R : struct with fields PX,PY,PZ,SX,SY,SZ (each nwindows x 1).
            obj.ensureObjectIsLoadedAndValid_();
            [sG, pG] = obj.build_geometry_libstruct_(G);
            [sE, pE] = obj.build_earth_libstruct_(E);
            [sR, pR] = obj.build_response_libstruct_();

            status = calllib(obj.LibAlias, 'forwardmodel', obj.HandleValue, pG, pE, pR);
            obj.throwIfNotOk_(status, 'forwardmodel');
            R = gatdaem1d_system_handle.unpack_response(sR);
        end
        function R = derivative(obj, G, E, dtype, dlayer)
            %DERIVATIVE Compute a derivative response as defined by the library's dtype strings.
            %
            %   R = S.derivative(G,E,dtype,dlayer)
            %
            % dtype is a string understood by the C++ library (see its possible-values message
            % when an invalid dtype is passed). Some derivative types require a layer index:
            %   - For conductivity derivatives (e.g. dtype='DC'), dlayer must satisfy
            %       0 <= dlayer <= nlayers-1
            %   - For thickness derivatives (e.g. dtype='DT'), dlayer must satisfy
            %       0 <= dlayer <= nlayers-2
            obj.ensureObjectIsLoadedAndValid_();
            [sG, pG] = obj.build_geometry_libstruct_(G);
            [sE, pE] = obj.build_earth_libstruct_(E);
            [sR, pR] = obj.build_response_libstruct_();

            dtypePtr = libpointer('cstring', char(dtype));
            dlayer   = int32(dlayer);

            status = calllib(obj.LibAlias, 'derivative', obj.HandleValue, pG, pE, dtypePtr, dlayer, pR);
            obj.throwIfNotOk_(status, 'derivative');

            R = gatdaem1d_system_handle.unpack_response(sR);
        end
        
        function A = fm_dlogc(obj, G, E)
            %FM_DLOGC Forward model plus d/d(log(conductivity)) for all layers in one call.
            %
            %   A = S.fm_dlogc(G,E)
            %
            % Output format
            %   A.FM        : forward response with fields PX,PY,PZ,SX,SY,SZ
            %   A.dlogC(k)  : derivative wrt log conductivity of layer k (1..nlayers),
            %                with the same fields as A.FM.
            % A = obj.fm_dlogc(G,E)
            %
            % Output buffer assumed packed as:
            %   A.FM:    [1 x 1 struct]       forward model - 6 blocks (PX,PY,PZ,SX,SY,SZ), each length nwindows
            %   A.dlogC: [nLayers x 1 struct] d/d(log-base-eC[k]) for each layer k=1..nlayers: in same channel order as FM
            
            obj.ensureObjectIsLoadedAndValid_();
            [sG, pG] = obj.build_geometry_libstruct_(G);
            [sE, pE] = obj.build_earth_libstruct_(E);

            nw = double(obj.nwindows());
            nl = double(sE.nlayers);

            nTotal = 6 * nw * (1 + nl);
            bufPtr = libpointer('doublePtr', zeros(nTotal,1));

            status = calllib(obj.LibAlias, 'fm_dlogc', obj.HandleValue, pG, pE, bufPtr);
            obj.throwIfNotOk_(status, 'fm_dlogc');

            v = bufPtr.Value;

            % Helper function: kth block (1-based) of length nw
            blk = @(k) v((k-1)*nw + (1:nw));

            % Forward model FM
            A.FM.PX = blk(1);
            A.FM.PY = blk(2);
            A.FM.PZ = blk(3);
            A.FM.SX = blk(4);
            A.FM.SY = blk(5);
            A.FM.SZ = blk(6);

            % Derivatives per layer
            A.dlogC = repmat(struct('PX',[],'PY',[],'PZ',[],'SX',[],'SY',[],'SZ',[]), nl, 1);
            for k = 1:nl
                b = 6 + (k-1)*6;
                A.dlogC(k).PX = blk(b+1);
                A.dlogC(k).PY = blk(b+2);
                A.dlogC(k).PZ = blk(b+3);
                A.dlogC(k).SX = blk(b+4);
                A.dlogC(k).SY = blk(b+5);
                A.dlogC(k).SZ = blk(b+6);
            end
        end
    end

    methods (Access = private)
        function ensureLibraryIsLoaded_(obj)
            if ~libisloaded(obj.LibAlias)
                %%This function is just about loading the library - its not the creation of the system handle
                [notfound, warnings] = loadlibrary(obj.DllPath, obj.ProtoFcn, 'alias', obj.LibAlias);
                if(~isempty(notfound))                  
                    error(warnings);
                end                
            end

            if isempty(obj.LibraryABIVersion)
                obj.LibraryABIVersion = calllib(obj.LibAlias, 'gatdaem1d_get_abi_version');                            
            end

            if isempty(obj.Enums)
                [~,~,enuminfo,~] = gatdaem1d_proto();
                obj.Enums.IPTYPE = enuminfo.GATDAEM1D_IPTYPE;
                obj.Enums.STATUS = enuminfo.GATDAEM1D_STATUS;
                obj.ProtoABIVersion = enuminfo.GATDAEM1D_EXPECTED_ABI.GATDAEM1D_EXPECTED_ABI_VERSION;                
                if(obj.LibraryABIVersion ~= obj.ProtoABIVersion)
                    msg = sprintf('The gatdaem1d library ABI version (%d) does not match the proto file ABI version (%d).\nPlease ensure that the loaded library (gatdaem1d.dll or .so) and the proto file (gatdaem1d_proto.m) are from the same release.',obj.LibraryABIVersion, obj.ProtoABIVersion);
                    error(msg);
                end
            end
        end
        function ensureObjectIsLoadedAndValid_(obj)
            ensureLibraryIsLoaded_(obj);
            if(~obj.IsValid)
                msg = sprintf("ensureObjectIsLoadedAndValid_() failed.\nDid you previously delete the system handle with free_system_handle()?\n%s",obj.SystemFile);
                error(msg);
            end            
        end
        function throwIfNotOk_(obj, status, fn)
            % GATDAEM1D_STATUS enums
            %              GATDAEM1D_STATUS_OK: 0
            %    GATDAEM1D_STATUS_ERRBADHANDLE: 1
            %     GATDAEM1D_STATUS_ERRBADMAGIC: 2
            %      GATDAEM1D_STATUS_ERRNULLSYS: 3
            % GATDAEM1D_STATUS_ERRCPPEXCEPTION: 4

            status = int32(status);            
            if status == obj.Enums.STATUS.GATDAEM1D_STATUS_OK
                return;
            end

            try
                cppmsg = gatdaem1d_get_last_error_string(obj, status);
            catch
                cppmsg = '';
            end

            if isempty(cppmsg)
                err_msg = sprintf('gatdaem1d call failed with status=%d in function %s, with message\n', status, fn);
            else
                err_msg = sprintf('gatdaem1d call failed with status=%d in function %s, with message\n%s\n', status, fn, cppmsg);
            end

            err_msg = sprintf('SystemFile is %s\n', obj.SystemFile);

            error(err_msg);
        end
        function [sR, pR] = build_response_libstruct_(obj)
            obj.ensureObjectIsLoadedAndValid_();
            nw = int32(obj.nwindows());
            sR = libstruct('gatdaem1d_response');
            sR.nwindows = nw;

            sR.PX = libpointer('doublePtr', zeros(double(nw),1));
            sR.PY = libpointer('doublePtr', zeros(double(nw),1));
            sR.PZ = libpointer('doublePtr', zeros(double(nw),1));
            sR.SX = libpointer('doublePtr', zeros(double(nw),1));
            sR.SY = libpointer('doublePtr', zeros(double(nw),1));
            sR.SZ = libpointer('doublePtr', zeros(double(nw),1));

            pR = libpointer('gatdaem1d_responsePtr', sR);
        end
        function [sG, pG] = build_geometry_libstruct_(obj,G)
            obj.ensureObjectIsLoadedAndValid_();
            sG = libstruct('gatdaem1d_system_geometry');
            f = fieldnames(sG);
            for i = 1:numel(f)
                if ~isfield(G, f{i})
                    error('build_geometry_libstruct_: missing field G.%s', f{i});
                end
                val = G.(f{i});
                if ~isscalar(val)
                    error('build_geometry_libstruct_: G.%s must be scalar', f{i});
                end
                if ~isreal(val)
                    error('build_geometry_libstruct_: G.%s must be real', f{i});
                end
                sG.(f{i}) = double(val);
            end
            pG = libpointer('gatdaem1d_system_geometryPtr', sG);
        end
        function [sE, pE] = build_earth_libstruct_(obj,E)
            % Matches structs.gatdaem1d_earth in proto:
            % iptype:int32, nlayers:int32, thickness:doublePtr, conductivity:doublePtr,
            % chargeability:doublePtr, timeconstant:doublePtr, frequencydependence:doublePtr
            obj.ensureObjectIsLoadedAndValid_();
            if ~isfield(E,'conductivity')
                error('build_earth_libstruct_: missing E.conductivity');
            end
            if ~isfield(E,'thickness')
                error('build_earth_libstruct_: missing E.thickness');
            end

            c = double(E.conductivity(:));
            t = double(E.thickness(:));

            nl = int32(numel(c));
            if numel(t) ~= max(0, double(nl)-1)
                error('build_earth_libstruct_: thickness must have length nlayers-1');
            end

            sE = libstruct('gatdaem1d_earth');
            sE.iptype       = int32(obj.Enums.IPTYPE.GATDAEM1D_IPTYPE_NONE);
            sE.nlayers      = int32(nl);
            sE.thickness    = libpointer('doublePtr', t);
            sE.conductivity = libpointer('doublePtr', c);

            % Optional IP arrays: pass NULL if absent/empty
            % GATDAEM1D_IPTYPE enums
            %     GATDAEM1D_IPTYPE_NONE: 0
            % GATDAEM1D_IPTYPE_COLECOLE: 1
            %   GATDAEM1D_IPTYPE_PELTON: 2
            hasIP = false;                    
            if(isfield(E,'iptype') && ~isempty(E.iptype))                
                if(E.iptype == obj.Enums.IPTYPE.GATDAEM1D_IPTYPE_COLECOLE || E.iptype == obj.Enums.IPTYPE.GATDAEM1D_IPTYPE_PELTON)
                    hasIP = true;                    
                    sE.iptype = int32(E.iptype);
                end
            end

            if(hasIP)                
                if isfield(E,'chargeability') && ~isempty(E.chargeability)
                    q = double(E.chargeability(:));
                    if numel(q) ~= double(nl), error('build_earth_libstruct_: chargeability must be length nlayers'); end
                    sE.chargeability = libpointer('doublePtr', q);
                else
                    error(['build_earth_libstruct_(): earth model has E.iptype as IPTYPE.GATDAEM1D_IPTYPE_COLECOLE or IPTYPE.GATDAEM1D_IPTYPE_PELTON but chargeability field is missing']);
                end

                if isfield(E,'timeconstant') && ~isempty(E.timeconstant)
                    tc = double(E.timeconstant(:));
                    if numel(tc) ~= double(nl), error('build_earth_libstruct_: timeconstant must be length nlayers'); end
                    sE.timeconstant = libpointer('doublePtr', tc);
                else
                    error(['build_earth_libstruct_(): earth model has E.iptype as IPTYPE.GATDAEM1D_IPTYPE_COLECOLE or IPTYPE.GATDAEM1D_IPTYPE_PELTON but timeconstant field is missing']);
                end

                if isfield(E,'frequencydependence') && ~isempty(E.frequencydependence)
                    fd = double(E.frequencydependence(:));
                    if numel(fd) ~= double(nl), error('build_earth_libstruct_: frequencydependence must be length nlayers'); end
                    sE.frequencydependence = libpointer('doublePtr', fd);
                else
                    error(['build_earth_libstruct_(): earth model has E.iptype as IPTYPE.GATDAEM1D_IPTYPE_COLECOLE or IPTYPE.GATDAEM1D_IPTYPE_PELTON but frequencydependence field is missing']);
                end
            else
                sE.chargeability       = libpointer('doublePtr', []); % true NULL
                sE.timeconstant        = libpointer('doublePtr', []); % true NULL
                sE.frequencydependence = libpointer('doublePtr', []); % true NULL
            end

            pE = libpointer('gatdaem1d_earthPtr', sE);
        end

    end %% methods (Access = private)

    methods (Static, Access = private)
        function R = unpack_response(sR)
            R.PX = sR.PX;
            R.PY = sR.PY;
            R.PZ = sR.PZ;
            R.SX = sR.SX;
            R.SY = sR.SY;
            R.SZ = sR.SZ;
        end
    end
end
