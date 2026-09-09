#!/usr/bin/env python3
"""Run with `sage -python OtherTools/export_xhash_test_suites.py`."""

import ast
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil

from sage.all import GF, Matrix, PolynomialRing, QQ, gcd, inverse_mod, vector

ROOT = Path(__file__).resolve().parents[1]
XHASH = ROOT / "Attack" / "Xhash"
P0 = 4611686018427388039
P1 = 4611686018427388073


def old_case(filename, name, width, free, degree, block, steps):
    """Read the literal original polynomials and triangular ideal without running it."""
    text = (XHASH / filename).read_text().replace("^", "**")
    tree = ast.parse(text)
    ring = PolynomialRing(QQ, names=[f"x{i}" for i in range(7)])
    env = dict(zip(ring.variable_names(), ring.gens()))
    polynomials = {}
    relations = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        if re.fullmatch(r"f\d+", target.id) or target.id == "ideal":
            value = eval(compile(ast.Expression(node.value), filename, "eval"), {"__builtins__": {}}, env)
            if target.id == "ideal":
                relations = value
            else:
                polynomials[int(target.id[1:])] = value
    equations = [polynomials[i] for i in range(len(polynomials))]
    return dict(name=name, p=P0, variables=[f"x{i}" for i in range(len(equations))],
                free_inputs=free, equations=list(map(str, equations)),
                relations=[[str(a), str(b)] for a, b in relations],
                parameters=dict(state_width=width, alpha=degree, nonlinear_block=block,
                                approximate_steps=steps, fixed_outputs=free,
                                input_state=[f"x{i}" for i in range(free)] + ["0"] * (width-free)),
                provenance=dict(source=filename, model="original XHash-style toy polynomial system",
                                zero_output_constraints=True, permutation_claim=False))


def unified_cases():
    field = GF(P1)
    assert gcd(3, P1 - 1) == gcd(3, P1**3 - 1) == 1
    univariate = PolynomialRing(field, "z")
    z = univariate.gen()
    modulus = next(z**3 + z + field(c) for c in range(1, 100)
                   if (z**3 + z + field(c)).is_irreducible())
    # Compute the coordinate polynomials in K[a,b,c][z]/(modulus).
    coordinates = PolynomialRing(field, names=["a", "b", "c"])
    a, b, c = coordinates.gens()
    extension_polys = PolynomialRing(coordinates, "z")
    zz = extension_polys.gen()
    lifted_modulus = sum(coordinates(modulus[i]) * zz**i for i in range(4))
    power = ((a + b*zz + c*zz**2)**3) % lifted_modulus
    coordinate_power = [power[i] for i in range(3)]
    mds = Matrix(field, [[1 / field(i + j + 4) for j in range(4)] for i in range(4)])
    for size in range(1, 5):
        for rows in itertools.combinations(range(4), size):
            for cols in itertools.combinations(range(4), size):
                assert mds.matrix_from_rows_and_columns(rows, cols).det() != 0
    seed = b"dixon-xhash-unified-v1-p4611686018427388073-t4-alpha3"
    raw = hashlib.shake_256(seed).digest(8 * 32)
    constants = [field(int.from_bytes(raw[8*i:8*i+8], "little")) for i in range(32)]
    inverse = int(inverse_mod(3, P1 - 1))
    metadata = dict(state_width=4, alpha=3, nonlinear_blocks=[3, 1],
                    extension_modulus=str(modulus),
                    coordinate_power=list(map(str, coordinate_power)),
                    matrix=[[int(v) for v in row] for row in mds.rows()],
                    constant_seed=seed.decode(), constants=list(map(int, constants)),
                    row_vector_convention=True, base_rate=2, base_capacity=2)
    result = []
    for name, free, steps in [("c1_3step_a3", 1, 3), ("c2_3step_a3", 2, 3), ("c1_5step_a3", 1, 5)]:
        count = free + (4 if steps == 3 else 8)
        ring = PolynomialRing(field, names=[f"x{i}" for i in range(count)])
        variables = ring.gens()
        equations, relations = [], []

        def mix(state):
            return list(vector(state) * mds)

        def fb_rhs(state, offset):
            state = mix([state[i] + constants[offset+i] for i in range(4)])
            state = mix([v**3 for v in state])
            return [state[i] + constants[offset+4+i] for i in range(4)]

        def p3f(state):
            state = [state[i] + constants[8+i] for i in range(4)]
            powered = [f(*state[:3]) for f in coordinate_power] + [state[3]**3]
            return mix([powered[i] + constants[12+i] for i in range(4)])

        initial = list(variables[:free]) + [ring(0)] * (4-free)
        first = list(variables[free:free+4])
        for lhs, rhs in zip([v**3 for v in first], fb_rhs(initial, 0)):
            relations.append((lhs, rhs))
            equations.append(lhs-rhs)
        middle = p3f(first)
        if steps == 3:
            equations.extend(middle[:free])
        else:
            second = list(variables[free+4:])
            for lhs, rhs in zip(second, middle):
                relations.append((lhs, rhs))
                equations.append(lhs-rhs)
            equations.extend(fb_rhs(second, 16)[:free])
        assert len(equations) == count
        # Independent finite-field extension oracle for the cubic coordinates.
        extension = field.extension(modulus, "u")
        u = extension.gen()
        for trial in range(5):
            numeric_input = [field(3 + trial*7 + i) for i in range(free)] + [field(0)] * (4-free)
            numeric_first = [v**inverse for v in fb_rhs(numeric_input, 0)]
            shifted = [numeric_first[i]+constants[8+i] for i in range(4)]
            numeric_extension = (shifted[0] + shifted[1]*u + shifted[2]*u**2)**3
            coefficients = list(numeric_extension.polynomial())
            coefficients += [field(0)] * (3-len(coefficients))
            reference = mix([coefficients[i]+constants[12+i] for i in range(3)]
                            + [shifted[3]**3+constants[15]])
            assert p3f(numeric_first) == reference
            values = numeric_input[:free] + numeric_first
            if steps == 5:
                values += reference
            assert all(f(*values) == 0 for f in equations[:len(relations)])
            expected = reference if steps == 3 else fb_rhs(reference, 16)
            assert [f(*values) for f in equations[len(relations):]] == expected[:free]
        result.append(dict(name=name, p=P1, variables=list(ring.variable_names()),
                           free_inputs=free, equations=list(map(str, equations)),
                           relations=[[str(a), str(b)] for a, b in relations],
                           parameters=dict(metadata, approximate_steps=steps, fixed_outputs=free,
                                           input_state=[str(v) for v in initial]),
                           provenance=dict(model="unified reduced toy permutation, not standard XHash12",
                                           zero_output_constraints=True, permutation_claim=True,
                                           validation="all Cauchy minors; power-map gcds; 5 forward checks per case")))
    return result


def export(directory, cases, description):
    from xhash_suite_runner import digest, read_case
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "data").mkdir(exist_ok=True)
    (directory / "_support").mkdir(exist_ok=True)
    shutil.copy2(ROOT / "OtherTools/xhash_suite_runner.py", directory / "_support/runner.py")
    shutil.copy2(ROOT / "drsolve/drsolve_sage_interface.sage", directory / "_support/drsolve_sage_interface.sage")
    manifest = []
    for data in cases:
        name = data["name"]
        data["sha256"] = digest(data)
        casepath = directory / "data" / (name + ".json")
        casepath.write_text(json.dumps(data, indent=2) + "\n")
        read_case(casepath)
        equations = data["equations"]
        variables = ",".join(data["variables"])
        magma = f'''// {name}; authoritative input SHA256: {data["sha256"]}
// Scope: grevlex Groebner basis only; excludes FGLM, roots and back-substitution.
SetSeed(20260908);
SetNthreads(1);
F := GF({data["p"]});
P<{variables}> := PolynomialRing(F, {len(data["variables"])}, "grevlex");
ps := [
    ''' + ",\n    ".join(equations) + f'''
];
printf "Case: {name}; equations: %o; input SHA256: {data["sha256"]}\\n", #ps;
t_cpu := Cputime();
t_wall := Realtime();
G := GroebnerBasis(ideal<P | ps>);
printf "GB CPU seconds: %o; wall seconds: %o; basis length: %o\\n", Cputime(t_cpu), Realtime(t_wall), #G;
'''
        (directory / (name + "_gb.magma")).write_text(magma)
        for method in ("direct", "hybrid"):
            script = f'''# Run: sage {name}_{method}.sage [--check-only] [--timeout 3600] [--threads 1]
from pathlib import Path
import sys
suite_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(suite_dir / "_support"))
from runner import run
run("{name}", "{method}", suite_dir)
'''
            (directory / (name + "_" + method + ".sage")).write_text(script)
        manifest.append({key: data[key] for key in ("name", "p", "parameters", "sha256")})
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    entries = "\n".join(f"| {d['name']} | {d['parameters']['approximate_steps']} | {d['parameters']['state_width']} | {d['free_inputs']} | {d['parameters']['alpha']} | {len(d['equations'])} |" for d in cases)
    first = cases[0]["name"]
    readme = f'''# {directory.name}

{description}

每个案例提供 3 个测试入口，共 9 个：`*_gb.magma`、`*_direct.sage`、`*_hybrid.sage`。
`data/` 是权威输入，`manifest.json` 记录参数，`_support/` 是公共运行依赖，不是额外测试。
同一案例的三个入口使用相同方程；SHA256 随输入和运行摘要保存。

| 案例 | 约步数 | 状态宽度 | 自由输入/输出约束数 | alpha | 方程/变量数 |
|---|---:|---:|---:|---:|---:|
{entries}

“约3步”对应 F,B,P3；“约5步”对应 F,B,P3,F,B。末端仿射处理不另计。
C1/C2 是本测试集的自由输入/输出约束数简称，不宣称与标准 CICO 定义无条件相同。
这些是缩小版或结构化 toy 实例，不是完整标准 XHash12。输出固定为零，未人为植入解；可能无解。

## 运行

从本目录运行（Sage 脚本也支持绝对路径启动）：

```bash
sage {first}_direct.sage --check-only
sage {first}_direct.sage --timeout 3600 --threads 1
sage {first}_hybrid.sage --timeout 3600 --threads 1
timeout 3600 magma {first}_gb.magma
```

其余案例替换文件名前缀即可。默认自动寻找仓库内的 `drsolve/drsolve`；也可加
`--drsolve /absolute/path/to/drsolve` 或设置 `DRSOLVE`。复制目录到其他机器时需要 Sage、
Magma，以及可运行的 DRSolve 二进制及其动态库；Sage 接口已随目录附带。
如果归档丢失二进制的执行权限，运行器会复制到本目录 `results/.bin/` 并为副本补权限。

## 计时范围和结果

- GB：grevlex 基计算，保留原实验的范围；不是包含 FGLM 的完整求解时间。
- direct：直接 Dixon 消去 x1,...，保留 x0。
- hybrid：对每个输出约束依次消去三角变量，使用 DixonIdeal 约化；C2 最后以两多项式 resultant 消去 x1。
- C1 也使用这条已有实现的约化路线，不将理论优化公式当作实测时间。
- Sage 默认单线程、总调用预算 3600 秒。输出保存在 `results/案例_方法_时间戳/`，
  包含逐阶段输入输出、`eliminant.txt` 和 `summary.json`。总时间包含调用开销和二进制内部可能执行的一元求根，
  不包含完整回代验证。失败/超时不会写成成功。
- Magma 输出到控制台；可以 `> case_gb.log 2>&1` 保存。它没有在当前导出环境中实际执行过。
- 不要把上述不同阶段的数字直接称为完整 CICO 求解时间。

重新生成：在仓库根目录运行 `sage -python OtherTools/export_xhash_test_suites.py`。
手动修改 JSON 后不要继续用旧 Magma 文件；应重新导出全部入口。
'''
    (directory / "README.md").write_text(readme)
    print(directory, "9 test entry points")


def main():
    minimal = [
        old_case("xhash-C1s3-m2.sage", "c1_5step_a3", 3, 1, 3, 3, 5),
        old_case("xhash-C1s3-m3.sage", "c1_5step_a7", 3, 1, 7, 3, 5),
        old_case("xhash-C2s2-m1+3.sage", "c2_3step_a3", 4, 2, 3, 2, 3),
    ]
    export(XHASH / "scheme1_minimal", minimal,
           "方案一：保留原方程与原域，修正文件配对。两个 C1 案例只改变指数；C2 保留原二维非线性块。原文件均未修改。")
    export(XHASH / "scheme2_unified", unified_cases(),
           "方案二：共同四维模板、alpha=3、域 p=4611686018427388073，真实三维扩域三次幂加一个基域三次幂。"
           "共同的 Cauchy MDS 与常量不随 C1/C2 改变。选择 C1/C2 约3步及 C1 约5步，共三个案例。")


if __name__ == "__main__":
    main()
