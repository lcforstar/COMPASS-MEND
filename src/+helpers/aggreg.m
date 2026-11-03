function [M] = aggreg(atmfile, SCOPEspec)

% 将 MODTRAN 大气辐射传输数据按 SCOPE 光谱波段进行矩形窗口平均聚合。
%
% 输入：
%   atmfile   - MODTRAN 输出的 .atm 文件路径（FLEX-S3_std.atm 格式）
%   SCOPEspec - SCOPE 光谱波段定义结构体（含 nreg, start, end, res）
%
% 输出：
%   M         - (nwS × 6) 矩阵，每列为对应 SCOPE 波段的平均大气参数：
%               列1: toasun  大气顶（TOA）太阳辐照度，用作 Eso*cos(tts)/pi [W/m²/μm]
%               列2: rdd     大气双向半球反射率（漫射→漫射），表示大气将
%                            向下漫射辐射重新反射回地面的比例 [无量纲 0-1]
%               列3: tss     太阳直射束透过率（beam→beam），直射光穿越大气
%                            到达地面且未发生散射的比例 [无量纲 0-1]
%               列4: tsd     太阳直射转漫射透过率（beam→diffuse），直射光
%                            被散射后以漫射形式到达地面的比例 [无量纲 0-1]
%               列5: tssrdd  预乘积 tss × rdd，用于修正地表与大气之间的
%                            多次散射耦合效应 [无量纲]
%               列6: Lab     大气背景下行辐亮度（sky background radiance），
%                            包含大气散射和热辐射贡献 [W/m²/sr/μm]
%
% .atm 文件列结构（T = 第3~20列，共18列）：
%   T列  文件变量    物理含义
%   ---  --------   --------
%    1   toasun      大气顶太阳辐照度
%    2   rso         大气路径反射率（太阳方向→传感器，后向散射贡献）
%    3   rdd         大气双向半球反射率（漫射↔漫射）
%    4   tss         直射太阳束透过率（beam→beam）
%    5   tsd         直射转漫射透过率（beam→diffuse）
%    6   too         上行直射透过率（地面→传感器，beam→beam）
%    7   tdo         上行漫射透过率（diffuse→传感器方向）
%    8   tsstoo      预乘积 tss × too
%    9   tsdtoo      预乘积 tsd × too
%   10   tsstdo      预乘积 tss × tdo
%   11   tsdtdo      预乘积 tsd × tdo
%   12   tssrdd      预乘积 tss × rdd（多次散射修正项）
%   13   toordd      预乘积 too × rdd
%   14   tssrddtoo   预乘积 tss × rdd × too
%   15   Lat         太阳方向大气程辐射（path radiance）
%   16   Lab         大气背景下行辐亮度
%   17   Labtoo      预乘积 Lab × too
%   18   Labtdo      预乘积 Lab × tdo
%
% 本函数仅提取 T 列 1,3,4,5,12,16 共 6 列（对应上表加粗项）。

% 读取 .atm 文件
s   = importdata(atmfile);
wlM = s.data(:, 2);          % 波长列 (nm)
T   = s.data(:, 3:20);       % 18列大气参数（T1~T18）

% 提取 6 个有效列：toasun(1), rdd(3), tss(4), tsd(5), tssrdd(12), Lab(16)
U     = [T(:, 1) T(:, 3) T(:, 4) T(:, 5) T(:, 12) T(:, 16)];

nwM   = length(wlM);

nreg  = SCOPEspec.nreg;    % 光谱区间数（SCOPE分为可见光、红外等多段）
streg = SCOPEspec.start;   % 各区间起始波长 (nm)
enreg = SCOPEspec.end;     % 各区间终止波长 (nm)
width = SCOPEspec.res;     % 各区间波段分辨率 (nm)

% 计算各区间的波段数
nwreg = int32((enreg - streg) ./ width) + 1;

% 各区间在输出矩阵中的起始偏移量
off   = int32(zeros(nreg, 1));

for i = 2:nreg
    off(i) = off(i - 1) + nwreg(i - 1);
end

nwS = sum(nwreg);
n   = zeros(nwS, 1);    % 每个 SCOPE 波段内累计的 MODTRAN 数据点数
S   = zeros(nwS, 6);    % 各波段 6 列参数的累加和

j   = int32(zeros(nreg, 1));  % 每个区间内当前波长对应的波段索引

for iwl = 1:nwM
    w   = wlM(iwl);    % 当前 MODTRAN 波长 (nm)
    for r = 1:nreg
        % 将 MODTRAN 波长映射到 SCOPE 第 r 区间的波段编号（1-based）
        j(r) = int32(round(w - streg(r)) ./ (width(r))) + 1;
        if j(r) > 0 && j(r) <= nwreg(r)      % 确认波长落在该区间范围内
            k       = j(r) + off(r);          % 对应的全局 SCOPE 波段索引
            S(k, :) = S(k, :) + U(iwl, :);   % 累加 6 列大气参数
            n(k)    = n(k) + 1;              % 计数加 1
        end
    end
end

% 各波段取平均，得到 SCOPE 光谱分辨率下的大气参数
M = zeros(size(S, 1), 6);
for i = 1:6
    M(:, i) = S(:, i) ./ n;
end

end
