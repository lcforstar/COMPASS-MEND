function fert_handle = initFertilizerController(fert_params)
%INITFERTILIZERCONTROLLER  Build fertilizer-event closure for N boundary injection.
%
% Returns a function handle fert_handle(t, layer) that returns NH4/NO3/urea
% mass inputs at a given time and layer when a fertilization event matches.
% Used by +init/setBoundaryCondition.m to wire fertilizer events into the
% N transport top boundary.
%
% Renamed from init_fertilizer_controller.m on 2026-05-13 as part of
% +soil/+nitrogen/ refactor (Skill §4 module organization).
%
% 输入 fert_params 字段:
%   t_fert      - datetime 向量，各施肥时间点
%   fert_kg_ha  - 1×n，各次施肥量 (kg/ha)
%   NH4_ratio   - 1×n，各次铵态氮比例 (0-1)
%   NO3_ratio   - 1×n，各次硝态氮比例 (0-1)
%   urea_ratio  - 1×n，各次尿素比例   (0-1)
%   depth_cm    - 标量，施肥层厚度 (cm)，= DEPTHEND - DEPTHSTART
%
% 返回句柄调用: [should_apply, F_conc, fR] = fert_handle(current_date)
%   should_apply - logical
%   F_conc       - 标量，当次总浓度
%   fR           - struct(.NH4_ratio, .NO3_ratio, .urea_ratio)，当次各形态比例

if nargin < 1, fert_params = struct(); end

n = numel(fert_params.t_fert);

% 默认值合并
if ~isfield(fert_params, 'NH4_ratio'),  fert_params.NH4_ratio  = zeros(1,n); end
if ~isfield(fert_params, 'NO3_ratio'),  fert_params.NO3_ratio  = zeros(1,n); end
if ~isfield(fert_params, 'urea_ratio'), fert_params.urea_ratio = ones(1,n);  end
if ~isfield(fert_params, 'depth_cm'),  fert_params.depth_cm  = 20;          end

% 广播标量比例到向量（fert_type 为单一类型时）
if isscalar(fert_params.NH4_ratio),  fert_params.NH4_ratio  = repmat(fert_params.NH4_ratio,  1, n); end
if isscalar(fert_params.NO3_ratio),  fert_params.NO3_ratio  = repmat(fert_params.NO3_ratio,  1, n); end
if isscalar(fert_params.urea_ratio), fert_params.urea_ratio = repmat(fert_params.urea_ratio, 1, n); end

state.t_fert      = fert_params.t_fert;
state.fert_kg_ha      = fert_params.fert_kg_ha;
state.NH4_ratio   = fert_params.NH4_ratio;
state.NO3_ratio   = fert_params.NO3_ratio;
state.urea_ratio  = fert_params.urea_ratio;
state.depth_cm    = fert_params.depth_cm;
state.applied     = false(1, n);
state.eps_tol     = hours(1);

fert_handle = @apply_fertilizer;

    function [should_apply, F_conc, fR] = apply_fertilizer(current_date)
        % Legacy stateful wrapper: the once-only flags live in this closure's
        % state.applied; the matching logic itself is the shared pure matcher
        % soil.nitrogen.fertilizer.matchFertilizerEvent (DRY). New-chain code calls
        % the pure matcher directly with runtime.management.fertilizer.applied.
        [should_apply, F_conc, fR, state.applied] = ...
            soil.nitrogen.fertilizer.matchFertilizerEvent(state, current_date, state.applied);
    end
end
