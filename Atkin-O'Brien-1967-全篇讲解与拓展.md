# 《p(n) 与 c(n) 模 13 的幂的同余性质》全篇讲解与拓展

**—— 面向高年级本科生的教学讲义**

> 论文：A. O. L. Atkin & J. N. O'Brien, *Some properties of p(n) and c(n) modulo powers of 13*, Trans. Amer. Math. Soc. **126** (1967), 442–459. 收稿日期 1966 年 5 月 3 日。

---

## 目录

1. [论文概览与核心结果](#1-论文概览与核心结果)
2. [基础知识补充](#2-基础知识补充)
3. [核心工具：模方程、π-赋值、U 算子](#3-核心工具)
4. [c(n) 部分：引理与定理的完整证明](#4-cn-部分)
5. [p(n) 部分：引理与定理的完整证明](#5-pn-部分)
6. [猜想 1、2 与其余定理](#6-猜想与其余定理)
7. [更高的视角：Hecke 理论、Serre 的模 p 模形式、现代发展](#7-更高的视角)
8. [附注与数值](#8-附注与数值)

---

## 1. 论文概览与核心结果

### 1.1 研究对象

论文研究两个经典数论序列模 $13$ 的幂的同余性质：

- **分拆函数** $p(n)$：把 $n$ 写成正整数之和的（无序）方法数，其生成函数为
  $$
  \frac{1}{f(x)}=\sum_{n=0}^{\infty}p(n)x^n,\qquad f(x)=\prod_{r=1}^{\infty}(1-x^r),\quad x=e^{2\pi i\tau}.
  $$

- **$j$-不变量的 Fourier 系数** $c(n)$：Klein 模不变量
  $$
  j(\tau)=\sum_{n=-1}^{\infty}c(n)x^n
  =\frac{\big(1+240\sum_{n\ge1}\sigma_3(n)x^n\big)^3}{x\,f^{24}(x)}-744,
  \qquad \sigma_3(n)=\sum_{d\mid n}d^3.
  $$
  展开式的前几项为
  $$
  j(\tau)=x^{-1}+744+196884x+21493760x^2+864299970x^3+\cdots,
  $$
  即 $c(-1)=1,\ c(0)=744,\ c(1)=196884,\ c(2)=21493760,\ c(3)=864299970,\dots$

### 1.2 问题的由来

Ramanujan 最早发现，对素数 $q\le 11$，存在一类"同余性质"：$p(5n+4)\equiv0\pmod 5$、$p(7n+5)\equiv0\pmod7$、$p(11n+6)\equiv0\pmod{11}$。Watson、Lehner、Atkin 等人将其推广为（论文中的 (1)、(2)）：

$$
(1)\qquad n\equiv0\pmod{2^a3^b5^c7^d11^e}\ \Longrightarrow\ c(n)\equiv0\pmod{2^{3a+8}3^{2b+5}5^c7^d11^e},
$$

$$
(2)\qquad 24n\equiv1\pmod{5^c7^d11^e}\ \Longrightarrow\ p(n)\equiv0\pmod{5^c7^{\lfloor(d+2)/2\rfloor}11^e}.
$$

**但对 $q>11$，这类"固定等差数列"的同余不再成立。** 于是 Newman 提出了两个更精细的问题（论文中的 (3)、(4)）：

> **(3)** 给定 $a,m$，$p(n)\equiv a\pmod m$ 是否对无穷多个 $n$ 可解？
>
> **(4)** 更强的：是否存在**正密度**的 $n$ 使 $p(n)\equiv a\pmod m$？

若能把 $p(n)$ 限制在某个**线性函数** $h(n)=bn+c$ 上（即 $p(bn+c)\equiv a\pmod m$ 对所有 $n$ 成立），则 (3) 与 (4) 同时成立。本文在 $m$ 为 $13$ 的幂时，对这两个问题做出了重要推进。

### 1.3 论文的主要结果

**定理 1**（$c(n)$ 的"递推"同余）。对每个 $a\ge1$，存在不被 $13$ 整除的常数 $k_a$，使对一切 $n$，
$$
c(13^{a+1}n)\equiv k_a\,c(13^a n)\pmod{13^a}.
$$
（$a=1$ 时即 Newman 的 $c(13^2n)\equiv8\,c(13n)\pmod{13}$，此时 $k_1\equiv8\pmod{13}$。）

**定理 2**（$p(n)$ 的"递推"同余）。定义
$$
P(N)=\begin{cases}p(n),&N=24n-1,\\0,&\text{其它情况}.\end{cases}
$$
则对每个 $a\ge1$，存在不被 $13$ 整除的常数 $K_a$，使对一切 $N$，
$$
P(13^{a+2}N)\equiv K_a\,P(13^a N)\pmod{13^a}.
$$

**定理 3**。$c(n)\equiv0\pmod{13^3}$ 对无穷多个 $n$ 成立。

**定理 4**。对每个 $a\ge1$ 与每个 $(a,13)=1$ 的 $a$，存在无穷多个 $n$ 使 $c(n)\equiv a\pmod{13^a}$。

**定理 7**。$p(n)\equiv0\pmod{13^4}$ 对无穷多个 $n$ 成立。

**定理 8**。对每个 $a\ge1$ 与每个 $(a,13)=1$ 的 $a$，存在无穷多个 $n$ 使 $p(n)\equiv a\pmod{13^a}$。

**定理 6**。$P(59^3\cdot13N)\equiv0\pmod{13}$ 当 $(N,59)=1$；等价地
$$
(14)\qquad p(59^4\cdot13n+111247)\equiv0\pmod{13}.
$$

**定理 9**。$P(13^2\cdot479\,n^2)\equiv0\pmod{13^2}$ 当 $(n,6)=1$；等价地
$$
p\!\Big(3373n^2-\frac{n^2-1}{24}\Big)\equiv0\pmod{13^2}.
$$

**定理 10**。$P(97^2\cdot103^2\cdot13^2N)\equiv0\pmod{13^2}$ 当 $\big(\frac N{97}\big)=\big(\frac N{103}\big)=-1$；从而
$$
(15)\qquad p(168544110546799\,n-6950975499605)\equiv0\pmod{13^2},\quad n\ge1.
$$

**定理 5**。猜想 2 在 $a=1,2$ 时成立。

**猜想 1**（Hecke 型乘法关系的提升）。设 $a\ge1$，令 $t(n)\equiv c(13^a n)/c(13^a)\pmod{13^a}$。则对素数 $p\neq13$，
$$
t(np)-t(n)t(p)+p^{11}t(n/p)\equiv0\pmod{13^a}.
$$

**猜想 2**（半整数权的 Hecke 型关系）。设 $a\ge1$，$p\neq13$ 为 $\ge5$ 的素数。则存在常数 $k=k(p,a)$，使对一切 $N$，
$$
P(p^2\cdot13^a N)-\Big\{k-\Big(\frac{-3N}{p}\Big)p^{12}\Big\}P(13^aN)+p^{-3}P(13^aN/p^2)\equiv0\pmod{13^a},
$$
其中 $\big(\tfrac{\cdot}{\cdot}\big)$ 为 Legendre 符号（$p^{-3}$ 指 $\bmod 13^a$ 的逆元）。

> **阅读提示**：定理 1、2 是"结构定理"，是所有其它结果的引擎。它们说的是：把 $c(n)$ 的下标乘以 $13$（或把 $P(N)$ 的下标乘以 $13^2$），相当于在值上乘以一个固定的 $13$-单位 $k_a$（或 $K_a$），且同余精度为 $13^a$。一旦有了这样的"递推同余"，再配合一个具体的"巧合零点"，就能产生无穷多个零点（定理 3、7），并能用群的循环结构填满所有单位剩余类（定理 4、8）。

---

## 2. 基础知识补充

### 2.1 整数分拆与 $p(n)$

**定义**。$n$ 的一个**分拆**是把 $n$ 写成一列非负不减正整数之和（次序无关）。$p(n)$ 表示 $n$ 的所有分拆数。约定 $p(0)=1$。

**生成函数**。每个分拆 $n=\lambda_1+\lambda_2+\cdots$（$\lambda_i$ 为部分，可重复）对应到因子展开
$$
\frac1{1-x}=\sum_{m\ge0}x^m
$$
中"每一部分 $m$ 可选任意重数"。于是
$$
\sum_{n\ge0}p(n)x^n=\prod_{m=1}^{\infty}\frac1{1-x^m}=\frac1{f(x)},\qquad f(x)=\prod_{m=1}^{\infty}(1-x^m).
$$

**Euler 五边形数定理**（贯穿全文的关键）：
$$
f(x)=\prod_{m=1}^\infty(1-x^m)=\sum_{k=-\infty}^{\infty}(-1)^k x^{k(3k-1)/2}.
$$
这给出 $p(n)$ 的一个（虽然不太快但完全初等的）递推公式，也是后文所有 $\eta$-商展开的基础。

**前几项**：
$$
p(0)=1,\ p(1)=1,\ p(2)=2,\ p(3)=3,\ p(4)=5,\ p(5)=7,\ p(6)=11,\ p(7)=15,\dots
$$
论文中会用到 $p(6)=11,\ p(19)=490,\ p(32)=8349,\ p(45)=89134$ 等值。

### 2.2 Dedekind $\eta$ 函数

设 $x=e^{2\pi i\tau}$，$\mathrm{Im}\,\tau>0$。Dedekind $\eta$ 函数定义为
$$
\eta(\tau)=e^{2\pi i\tau/24}\prod_{m=1}^{\infty}(1-e^{2\pi i m\tau})=x^{1/24}f(x).
$$
它的两个关键性质：

1. **模变换**：$\eta(-1/\tau)=\sqrt{-i\tau}\,\eta(\tau)$；更一般地，$\eta$ 是权 $\tfrac12$、带乘子系的模形式。这说明 $\eta$ 的商是模函数。
2. **模方程**：对素数 $\ell$，$\eta(\ell\tau)/\eta(\tau)$ 与 $\eta(\tau)$ 之间满足一个次数 $\ell$ 的代数方程（见 §3）。

论文里反复出现的是下面两个商：
$$
\varphi(x)=\frac{\eta(169\tau)}{\eta(\tau)}=x^7\,\frac{f(x^{169})}{f(x)},
\qquad
\psi(x)=\frac{\eta^2(13\tau)}{\eta^2(\tau)}=x\,\frac{f^2(x^{13})}{f^2(x)}.
\tag{16}
$$
（验证：$\eta(169\tau)/\eta(\tau)=x^{169/24}f(x^{169})/(x^{1/24}f(x))=x^{168/24}\frac{f(x^{169})}{f(x)}=x^7\frac{f(x^{169})}{f(x)}$；$\eta^2(13\tau)/\eta^2(\tau)=x^{26/24}f^2(x^{13})/(x^{2/24}f^2(x))=x\frac{f^2(x^{13})}{f^2(x)}$。）

### 2.3 模形式入门（只需最小限度的知识）

**定义（非严格但够用）**。一个（整权 $k$、级 $1$ 的）**模形式**是上半平面上的全纯函数 $F(\tau)$，满足
$$
F\!\Big(\frac{a\tau+b}{c\tau+d}\Big)=(c\tau+d)^k F(\tau),\qquad \begin{pmatrix}a&b\\c&d\end{pmatrix}\in SL_2(\mathbb Z),
$$
且在 $\infty$ 处全纯（即 $q$ 展开 $F=\sum a(n)q^n$ 无非正指数项）。若 $a(0)=0$ 则称**尖点形式**。

**两个原型**（$q=x$）：

- **Eisenstein 级数**（权 4）：$E_4(\tau)=1+240\sum_{n\ge1}\sigma_3(n)q^n$。
- **Ramanujan $\Delta$**（权 12 的唯一尖点形式，差常数倍）：
$$
\Delta(\tau)=\eta^{24}(\tau)=x\prod_{m\ge1}(1-x^m)^{24}=x f^{24}(x)=\sum_{n\ge1}\tau(n)x^n.
$$
$\tau(n)$ 是 **Ramanujan $\tau$ 函数**：$\tau(1)=1,\ \tau(2)=-24,\ \tau(3)=252,\dots$

**$j$-不变量**（权 $0$ 的模函数，即 $SL_2(\mathbb Z)$ 不变的亚纯函数）：
$$
j(\tau)=\frac{E_4^3}{\Delta}=\frac{(1+240\sum\sigma_3(n)q^n)^3}{q\prod(1-q^m)^{24}}.
$$
在 $\infty$ 处有单极点（$j=q^{-1}+744+O(q)$），故 $j-744$ 是"权 $0$、只在 $\infty$ 有极点"的模函数。其 Fourier 系数就是 $c(n)$。

**Hecke 特征值（关键事实）**。$\Delta=\sum\tau(n)q^n$ 是权 $12$ 的 Hecke 特征形式，故对素数 $p$：
$$
\tau(pn)=\tau(p)\tau(n)-p^{11}\tau(n/p).
\tag{$\star$}
$$
（$n/p$ 非整数时 $\tau(n/p)=0$。）这正是后文猜想 1 的"原型"。

### 2.4 Ramanujan 同余式的背景

Ramanujan（1919）发现并猜想，Watson（1938）、Atkin 等人证明：
$$
p(5n+4)\equiv0\pmod5,\quad p(7n+5)\equiv0\pmod7,\quad p(11n+6)\equiv0\pmod{11},
$$
且可推广到更高的幂（本文 (2)）。这些同余的深刻之处在于：**它们来自模形式的"整体"结构（Hecke 算子、Galois 表示），而非 $p(n)$ 的数值巧合**。本文正是用 $\eta$ 的模方程把这一思想推进到 $q=13$。

### 2.5 Newman 的两个问题

Newman（1960）提出 (3)：给定 $a,m$，是否总有 $p(n)\equiv a\pmod m$ 无穷多次？并进一步问 (4)：是否正密度。本文在 $m=13^a$ 时：

- 定理 4、8 肯定地回答了 (3)（对所有与 $13$ 互素的剩余类 $a$）。
- (4)（正密度）在本文中**仍未解决**（见 §7 的现代发展）。

### 2.6 13-adic 赋值 $\pi$

全文的核心记法。对整数 $a$，定义 $\pi(a)$ 为满足 $13^{\pi(a)}\Vert a$ 的指数（即 $13^{\pi(a)}$ 恰好整除 $a$），并约定 $\pi(0)=\infty$。对有理数 $a=b/c$，定义 $\pi(a)=\pi(b)-\pi(c)$。则（论文 (19)、(20)）：

$$
(19)\ \ \pi(ab)=\pi(a)+\pi(b),\qquad
(20)\ \ \pi(a+b)\ge\min\{\pi(a),\pi(b)\},
$$

且当 $\pi(a)\ne\pi(b)$ 时 (20) 取等号（这是"强三角不等式"，即 **13-adic 赋值** $v_{13}$ 的标准性质）。

**直觉**：$\pi(a)$ 度量 $a$ 被 $13$ 整除的"深度"。例如 $\pi(169)=2$、$\pi(5)=0$。全文所有引理都是在给某些系数 $a_{rp}$ 的 $\pi$ 值**下界**，从而得到"这些系数被 $13$ 的幂整除"的同余。

---

## 3. 核心工具

### 3.1 模方程（O'Brien 的结果 [1]）

**定理（模方程）**。设 $t=\varphi(x^{1/13})$、$g=\psi(x)$。则 $t,g$ 满足
$$
t^{13}+\beta_1 t^{12}+\cdots+\beta_{12}t+\beta_{13}=0,
\tag{17}
$$
其中
$$
\beta_r=\sum_{a}\beta_{ra}\,g^a,\qquad \beta_{ra}\in\mathbb Z,
$$
且（**Lemma 1**，由附录 C 的表直接验证）
$$
\pi(\beta_{ra})\ge\Big\lfloor\frac{13a-7r+13}{14}\Big\rfloor.
\tag{L1}
$$

**根的结构**。把 (17) 看作关于 $t$ 的方程（$g$ 固定），其 13 个根为
$$
t=\varphi(\omega^m x^{1/13}),\qquad m=0,1,\dots,12,
$$
其中 $\omega$ 是 13 次本原单位根。这是 (17) 成立的"理由"：方程 $t^{13}+\dots=0$ 的根正是 $\varphi$ 在 $x^{1/13}$ 的 13 个分歧点上的取值。

**为什么是 13 次、为什么是这些权？** 因为 $\varphi(x^{1/13})\sim x^{7/13}$（主导项），而 $\psi\sim x$，故单项式 $g^a t^{13-r}\sim x^{a+7(13-r)/13}$。方程齐次性要求 $a+7(13-r)/13=7$，即 $13a=7r$。所以 $\beta_r$ 中 $g$ 的幂次 $a$ 围绕 $7r/13$ 取值，这也解释了 Lemma 1 中"$13a-7r$"的出现。

### 3.2 幂和 $S_r$ 与 Newton 恒等式

记 $S_r$ 为 (17) 的 13 个根的 $r$ 次幂和。由根的结构与单位根求和（**Lemma 3**，见下），
$$
S_r=\sum_{m=0}^{12}\varphi^r(\omega^m x^{1/13})=13\,U\,\varphi^r(x),
\tag{29}
$$
即
$$
U\varphi^r(x)=\tfrac1{13}S_r.
$$

另一方面，$S_r$ 作为 $\psi$ 的多项式，可写
$$
S_r=\sum_p a_{rp}\,\psi^p,
$$
其系数 $a_{rp}$ 由 **Newton 恒等式** 从 $\beta_{ra}$ 递归决定。Newton 恒等式（对多项式 $t^{13}+\beta_1t^{12}+\cdots+\beta_{13}$）给出
$$
S_R=-\sum_{r=1}^{R-1}\beta_{R-r}S_r-R\beta_R\qquad(1\le R\le13).
\tag{$\dagger$}
$$

**Lemma 2**。$\pi(a_{rp})\ge\big\lfloor\frac{13p-7r+13}{14}\big\rfloor$，且 $a_{rp}=0$ 除非 $\lfloor(7r+12)/13\rfloor\le p\le 7r$。

**证明思路（全文的归纳模板）**。$R=1$ 时由 Lemma 1 直接得到（$S_1=-\beta_1$）。设对 $r<R$ 成立。由 ($\dagger$)，$a_{Rp}$ 是形如 $a_{r\sigma}\beta_{R-r,\,p-\sigma}$ 的整系数线性组合，故
$$
\pi(a_{Rp})\ge\min_{r,\sigma}\Big\{\pi(a_{r\sigma})+\pi(\beta_{R-r,p-\sigma})\Big\}.
$$
代入两个归纳假设，并用基本不等式（论文 (24)）
$$
\Big\lfloor\frac{A}{14}\Big\rfloor+\Big\lfloor\frac{B}{14}\Big\rfloor\ge\Big\lfloor\frac{A+B-13}{14}\Big\rfloor,
$$
即得
$$
\pi(a_{Rp})\ge\Big\lfloor\frac{13p-7R+13}{14}\Big\rfloor.
$$
归纳完成。**这是全文反复使用的"取 $\min$、化 $\lfloor\cdot\rfloor$ 求和"的技巧**，务必掌握。

### 3.3 U 算子

对任意幂级数 $F(x)=\sum A(n)x^n$，定义
$$
UF(x)=\sum_{n}A(13n)\,x^n.
\tag{27}
$$
它是**线性**的，并满足（论文 (28)，直接验证）：
$$
U\{F_1(x^{13})F_2(x)\}=F_1(x)\cdot UF_2(x).
\tag{28}
$$
（把 $F_1$ 的 $x^{13}$ 换成 $x$，把 $F_2$ 换成 $UF_2$。）

**Lemma 3（单位根筛法）**。设 $\omega\ne1,\ \omega^{13}=1$。则
$$
13\,UF(x)=\sum_{m=0}^{12}F(\omega^m x^{1/13}).
$$
（证明：右端 $=\sum_n A(n)x^{n/13}\sum_m\omega^{mn}$，而 $\sum_m\omega^{mn}=13$ 当 $13\mid n$，否则 $=0$。）

**推论**：$U\varphi^r=\frac1{13}S_r$，即 (29)。

**Lemma 4（U 在 $\psi^k$ 与 $\varphi\psi^k$ 上的作用）**。对 $k\ge1$，
$$
U\psi^k(x)=\sum_{r\ge1}C_{kr}\psi^r(x),\qquad
U\{\varphi(x)\psi^k(x)\}=\sum_{r\ge1}d_{kr}\psi^r(x),
\tag{30,31}
$$
其中 $C_{kr},d_{kr}\in\mathbb Z$ 满足
$$
\pi(C_{kr})\ge\Big\lfloor\frac{13r-k-1}{14}\Big\rfloor,\qquad
\pi(d_{kr})\ge\Big\lfloor\frac{13r-k-8}{14}\Big\rfloor,
\tag{32,33}
$$
且 $\pi(C_{11})=\pi(d_{11})=0$；$C_{kr}=0$ 除非 $\lfloor(k+12)/13\rfloor\le r\le13k$，$d_{kr}=0$ 除非 $\lfloor(k+19)/13\rfloor\le r\le13k+7$。

**为什么成立**：由 (29) 与 $S_r$ 的结构。例如
$$
U\psi^k=U\big\{x^{13k}f^{2k}(x^{13})\cdot f^{-2k}(x)\big\}
$$
经 (28) 与 $\psi$、$\varphi$ 的关系化为 $S_{2k}$ 的线性组合，再按 $S_r=\sum a_{rp}\psi^p$ 展开，即得 $C_{kr}$ 的 $\pi$-界由 Lemma 2 推出。

> **要点**：Lemma 4 说 $U$ 把"$\psi$ 的语言"映回"$\psi$ 的语言"，且每一步都"抬升" $13$-深度。这是下面两个归纳（$c(n)$ 与 $p(n)$）能够进行的原因。

---

## 4. c(n) 部分

### 4.1 Lemma 5（出发点）

$$
\sum_{n\ge1}c(13n)x^n=-\psi(x)+13^2\,U^2\psi(x).
\tag{L5}
$$

**证明**（附录 A，唯一用模理论的地方）。把 $g=\psi$ 看作 $\tau$ 的函数 $g(\tau)$，则 $g(\tau)$ 是 $\Gamma_0(13)$ 上的整函数（无极点）。由 Atkin [5] 的引理，
$$
13\{g(-1/13\tau)+13Ug(\tau)\}=g^{-1}(\tau)+13^2Ug(\tau)
$$
是 $\Gamma(1)$ 上的整函数。而 $g^{-1}$ 在 $i\infty$ 有单极点（留数 1），故 $j(\tau)=g^{-1}(\tau)+13^2Ug(\tau)+$ 常数。再证 $g^{-1}(-1/13\tau)+13Ug^{-1}(\tau)$ 是 $\Gamma(1)$ 上整函数且正则，从而是常数，得 $Ug^{-1}(\tau)=-g(\tau)+$ 常数。于是
$$
Uj(\tau)=-g(\tau)+13^2U^2g(\tau)+常数.
$$
由于 $Uj(x)=\sum c(13n)x^n$，即得 (L5)。$\blacksquare$

**推论**（与 $\tau(n)$ 的联系）。由 (L5) 模 13：
$$
c(13n)\equiv-[x^n]\psi(x)\pmod{13}.
$$
而 $f(x^{13})\equiv f(x)^{13}\pmod{13}$（Frobenius/"新手梦"），故
$$
\psi(x)=x\frac{f^2(x^{13})}{f^2(x)}\equiv x\frac{f^{26}(x)}{f^2(x)}=xf^{24}(x)=\Delta(x)=\sum\tau(n)x^n\pmod{13}.
$$
于是 $c(13n)\equiv-\tau(n)\pmod{13}$，这正是 Newman 证明 (6)、(7) 的出发点。

### 4.2 Lemma 6（$J_a$ 的 $\psi$ 展开）

记 $J_a(x)=\sum_{n\ge1}c(13^a n)x^n$。反复应用 Lemma 4 得
$$
J_a(x)=\sum_{r\ge1}j_{ar}\,\psi^r(x)\qquad(a\ge1),
\tag{35}
$$
其中 $r$ 从 $1$ 跑到 $13^{a+1}$，且
$$
\pi(j_{ar})\ge\Big\lfloor\frac{13r-2}{14}\Big\rfloor,\qquad \pi(j_{a1})=0.
\tag{36,37}
$$

**证明**。因为 $J_{a+1}=U J_a$（$U$ 的作用正是"取 $13$ 的倍数项"），故 $j_{a+1,r}=\sum_p j_{ap}C_{pr}$，由 (28)+(30) 得 $\pi(j_{a+1,r})\ge\min_p\{\pi(j_{ap})+\pi(C_{pr})\}$。代入归纳假设 (36) 与 (32)，取 $\min$ 在 $p=1$ 或 $2$ 处达到，即得 (36) 对 $a+1$。$\pi(j_{a1})=0$ 的保持类似。基情形 $a=1$ 由 Lemma 5 给出。$\blacksquare$

### 4.3 Lemma 7 与 Theorem 1

定义"交差行列式"（2×2 行列式）
$$
Y_{rs}^a=j_{a+1,r}\,j_{as}-j_{ar}\,j_{a+1,s}.
\tag{40}
$$

**Lemma 7**。$\pi(Y_{rs}^a)\ge a+\big\lfloor\frac{13(r+s)-32}{14}\big\rfloor$；特别地
$$
\pi(Y_{rs}^a)>a.
\tag{43}
$$

（证明仍是"取 $\min$"归纳：$Y^{a+1}=\sum_{p,\sigma}Y_{p\sigma}^a C_{pr}C_{\sigma s}$，用 (32) 与 (24) 式，$\min$ 在 $p+\sigma=3$ 达到。）

**Theorem 1 的证明**。由 (37)，$j_{a1}$ 是 $13$-单位，故可选 $k_a$（$13\nmid k_a$）使
$$
j_{a+1,1}\equiv k_a\,j_{a1}\pmod{13^a}.
$$
对 (43) 取 $s=1$：$j_{a+1,r}j_{a1}-j_{ar}j_{a+1,1}\equiv0\pmod{13^{a+1}}$。除以单位 $j_{a1}$，得
$$
j_{a+1,r}\equiv j_{ar}\,\frac{j_{a+1,1}}{j_{a1}}\equiv k_a\,j_{ar}\pmod{13^a}\quad(\forall r).
$$
故 $J_{a+1}(x)\equiv k_aJ_a(x)\pmod{13^a}$，比较 $x^n$ 系数即得
$$
c(13^{a+1}n)\equiv k_a\,c(13^a n)\pmod{13^a}.\qquad\blacksquare
$$

**最优性**：模数 $13^a$ 不能再改进为 $13^{a+1}$（由 $\pi(Y^a_{21})=a$ 看出，见论文 (44)）。

### 4.4 Theorem 3 与 Theorem 4

**Theorem 3**：论文发现（计算巧合）$c(5299\cdot13^3)\equiv0\pmod{13^3}$（其中 $5299=7\cdot757$）。由 Theorem 1 递推：
$$
c(5299\cdot13^{m+3})\equiv k_{m+3}\cdots k_4\,c(5299\cdot13^3)\equiv0\pmod{13^3}\quad(m\ge0).
$$
这些 $n=5299\cdot13^{m+3}$ 两两不同，故 $c(n)\equiv0\pmod{13^3}$ 无穷多次。$\blacksquare$

**Theorem 4**：需要理解 $k_a$ 的乘法群结构。计算得 $k_a\equiv k_2\pmod{13^2}$（$a\ge2$），且 $k_2\equiv41^3\pmod{13^2}$，而 $41$ 是 $13^2$ 的原根。因此 $k_a$ 在乘法群 $(\mathbb Z/13^a\mathbb Z)^\times$（阶 $12\cdot13^{a-1}$）中生成一个**指标 3 的子群**（阶 $4\cdot13^{a-1}$）。于是由 Theorem 1，若 $c(13^a n_0)$ 落在某个陪集中，则
$$
c(13^{a+m}n_0)\equiv k_a^m\,c(13^a n_0)\pmod{13^a}
$$
填满该陪集无穷多次。而直接计算给出三个不同的陪集代表：
$$
c(13\cdot1)\equiv-1,\quad c(13\cdot2)\equiv-2,\quad c(13\cdot5)\equiv6\pmod{13},
$$
这三个值落在指标 3 子群的三个不同陪集中，故并起来覆盖所有 $13$-单位剩余类。$\blacksquare$

> **核心思想**：Theorem 1 把"乘以 13 的幂"转化为"乘以 $k_a$"，而 $k_a$ 的循环作用 + 群结构 = 填满陪集。这是"递推同余 $\Rightarrow$ 分布结果"的经典模式。

---

## 5. p(n) 部分

### 5.1 记号 $L$ 与 Lemma 8

定义（论文 (46)，此处按正确形式书写）：
$$
L_{2a-1}(x)=f(x^{13})\sum_{n\ge1}P\big(13^{2a-1}(24n-13)\big)x^n,
$$
$$
L_{2a}(x)=f(x)\sum_{n\ge1}P\big(13^{2a}(24n-1)\big)x^n.
$$

**Lemma 8**：
$$
(47)\ L_1(x)=U\varphi(x),\qquad
(48)\ L_{2a}(x)=UL_{2a-1}(x),\qquad
(49)\ L_{2a+1}(x)=U\big(\varphi(x)L_{2a}(x)\big).
$$

**验证 (47)**（这是理解一切的钥匙）。直接计算：
$$
U\varphi=U\big[x^7f(x^{169})/f(x)\big]=f(x^{13})\,U\big[x^7/f(x)\big]
$$
（用 (28)，因 $f(x^{169})=f((x^{13})^{13})$ 是 $x^{13}$ 的函数）。而
$$
U\big[x^7/f(x)\big]=\sum_n p(13n-7)x^n,
$$
故 $U\varphi=f(x^{13})\sum p(13n-7)x^n$。又 $p(13n-7)=P(13(24n-13))$，得 $L_1=U\varphi$。$\blacksquare$

（(48)、(49) 类似地用 (28) 验证，(48) 是"即得"的。）

### 5.2 Lemma 9、10 与 Theorem 2

反复应用 Lemma 4，把 $L$ 写成 $\psi$ 的展开（论文 (50)、(51)）：
$$
L_{2a-1}(x)=\sum_r k_{ar}\psi^r(x),\qquad L_{2a}(x)=\sum_r l_{ar}\psi^r(x),
$$
其中（**Lemma 9**）
$$
\pi(k_{ar})\ge\Big\lfloor\frac{13r-9}{14}\Big\rfloor,\qquad
\pi(l_{ar})\ge\Big\lfloor\frac{13r-2}{14}\Big\rfloor,\qquad
\pi(k_{a1})=\pi(l_{a1})=0.
\tag{54,55,56}
$$

再定义 2×2 行列式（论文 (57)）
$$
\theta_{rs}^a=k_{a+1,r}k_{as}-k_{ar}k_{a+1,s},\qquad
\epsilon_{rs}^a=l_{a+1,r}l_{as}-l_{ar}l_{a+1,s}.
$$

**Lemma 10**：
$$
\pi(\theta_{rs}^a)>2a-1,\qquad \pi(\epsilon_{rs}^a)>2a.
\tag{64,65}
$$

**Theorem 2 的证明**。由 (56)，$k_{a1}$ 是单位，取 $K_{2a-1}$（$13\nmid K_{2a-1}$）使
$$
k_{a+1,1}\equiv K_{2a-1}\,k_{a1}\pmod{13^{2a-1}}.
$$
由 (64)（取 $s=1$）得 $k_{a+1,r}\equiv K_{2a-1}k_{ar}\pmod{13^{2a-1}}$，故
$$
L_{2a+1}(x)\equiv K_{2a-1}L_{2a-1}(x)\pmod{13^{2a-1}}.
$$
两边除以 $f(x^{13})$（合法，因 $f(x^{13})$ 的常数项为 1）并比较系数，得
$$
P(13^{2a+1}N)\equiv K_{2a-1}P(13^{2a-1}N)\pmod{13^{2a-1}}.
$$
（偶 $a$ 情形用 $\epsilon$ 完全类似。）令 $K_a=K_{2a-1}$（奇）或 $K_{2a}$（偶），即得 Theorem 2。$\blacksquare$

### 5.3 Theorem 7 与 Theorem 8

**Theorem 7**：计算巧合 $P(13^4\cdot22655)\equiv0\pmod{13^4}$。由 Theorem 2 递推
$$
P(13^{4+2m}\cdot22655)\equiv0\pmod{13^4}\quad(m\ge0),
$$
故 $p(n)\equiv0\pmod{13^4}$ 无穷多次。$\blacksquare$

**Theorem 8**：$K_a\equiv K_2\pmod{13^2}$（$a\ge2$），而 $K_2=45$ 是 $13^2$ 的原根，故 $K_a$ 是 $13^a$ 的原根（阶 $12\cdot13^{a-1}$）。因 $k_{a1},l_{a1}$ 是单位，$P(13^{2a-1}\cdot11)$、$P(13^{2a}\cdot23)$ 不被 $13$ 整除。于是
$$
P(13^{2a-1+2m}\cdot11)\ (\bmod 13^{2a-1})\quad\text{填满所有 }13\text{-单位剩余类},
$$
$$
P(13^{2a+2m}\cdot23)\ (\bmod 13^{2a})\quad\text{填满所有 }13\text{-单位剩余类}.
$$
故定理 8 成立。$\blacksquare$

### 5.4 Newman 序列（论文末尾的观察）

定义 $t(m)=P(13^{2m-1}N_0)$。Theorem 2 说明
$$
\pi(t(m+1))=\pi(t(m))\ \text{若}\ \pi(t(m))<2m-1,\quad
\pi(t(m+1))\ge2m-1\ \text{否则}.
$$
于是要么 $\pi(t(m))\ge2m-1$ 对所有 $m$（即"越来越深"），要么 $\pi(t(m))$ 最终稳定在某个 $<2m_0-1$ 的值。作者认为后者更可能——即这些序列的 $13$-深度通常**有界**。

---

## 6. 猜想与其余定理

### 6.1 猜想 1 与 Hecke 理论

猜想 1 说：规范化后的系数 $t(n)=c(13^a n)/c(13^a)\pmod{13^a}$ 满足
$$
t(np)-t(n)t(p)+p^{11}t(n/p)\equiv0\pmod{13^a}.
$$
**这就是权 12 Hecke 特征形式的乘法关系 ($\star$)**，只不过把 $\tau(n)$ 换成了 $c(13n)$（$\equiv-\tau(n)\bmod13$）再逐级提升。作者的观点是：Newman 的 (6)、(7) **本质上是 $\tau(n)$ 的 Hecke 乘法性在 $c(13n)$ 上的反映**，而非偶然。猜想 1 断言这一乘法性可以提升到模 $13^a$ 的每个精度。

### 6.2 猜想 2 与半整数权

$p(n)$ 对应的 $\eta^{-1}$ 是**半整数权**（权 $-\tfrac12$），其 Hecke 理论带有 Legendre 符号（二次扭）。猜想 2 是相应 Hecke 关系的模 $13^a$ 版本：
$$
P(p^2\cdot13^a N)-\Big\{k(p,a)-\Big(\frac{-3N}{p}\Big)p^{12}\Big\}P(13^aN)+p^{-3}P(13^aN/p^2)\equiv0\pmod{13^a}.
$$
其中 $\big(\frac{-3N}{p}\big)$ 的出现正是半整数权形式的特征（与 $\eta$ 的乘子系、与 $\theta$-对应相关）。**定理 5**：当 $a=1,2$ 时猜想 2 成立，其证明（§5）把问题化归到 Newman 关于 $P_{11},P_{23}$ 的已知结果 (68)、(69)，并用到 $P_{11}(13N)\equiv8P_{23}(N)\pmod{13}$ 等具体同余。

### 6.3 定理 6、9、10（由"巧合"与乘法性构造显式同余）

- **定理 6**：$P(59^3\cdot13N)\equiv0\pmod{13}$，用到了 $k(p,1)\equiv0\pmod{13}$ 对 $p=59$ 这一"巧合"。
- **定理 9**：$P(13^2\cdot479n^2)\equiv0\pmod{13^2}$（$(n,6)=1$），由 Lemma 11（$N_0=479$ 平方自由、$P(13^2N_0)\equiv0\pmod{13^2}$ 推出 $P(13^2n^2N_0)\equiv0\pmod{13^2}$）。
- **定理 10**：取最初两个使 $p^{2k}\equiv-1\pmod{13^2}$ 的素数 $97,103$，由猜想 2（$a=2$ 已证）推出 $P(97^2\cdot103^2\cdot13^2N)\equiv0\pmod{13^2}$ 当 $\big(\frac N{97}\big)=\big(\frac N{103}\big)=-1$，从而得到线性同余 (15)。

---

## 7. 更高的视角

### 7.1 Hecke 算子与 U 算子的统一理解

本文的 $U=U_{13}$ 算子与 Hecke 算子 $T_{13}$ 密切相关：对权 $k$ 形式 $F=\sum a(n)q^n$，
$$
T_{13}F=\sum_n\big(a(13n)+13^{k-1}a(n/13)\big)q^n,\qquad UF=\sum_na(13n)q^n.
$$
**关键事实**：$\Delta$ 是 $T_p$ 的特征形式（特征值 $\tau(p)$），故 $T_{13}\Delta=\tau(13)\Delta$。而本文的 Lemma 5 本质上是 $j$ 在 $U$ 作用下的分解。$U$ 作用在 $\eta$-商生成的子空间上是"有限维"的（Lemma 4 的 $C_{kr},d_{kr}$ 是有限个），这使归纳得以进行。

### 7.2 Serre 的模 $p$ 模形式理论与"乘法性"

猜想 1 的一个**现代解释**：模 $13$ 的模形式空间 $M_k(\mathbb F_{13})$ 中，$\bar\Delta=\Delta\bmod13$ 是 $T_p$ 的特征形式。Serre（1973）与 Swinnerton-Dyer（1973）建立了模 $\ell$ 模形式的系统理论，说明 $j$ 的系数 $c(n)\bmod\ell$ 受 $\Delta$ 的 Galois 表示控制。本文的猜想 1 可以理解为：$c(13^a n)/c(13^a)$ 在模 $13^a$ 意义下"继承"了 $\bar\Delta$ 的 Hecke 特征值结构。这正是后来 **Galois 表示** 与 **$p$-adic 模形式** 理论的早期萌芽。

### 7.3 与 Galois 表示的联系（简述）

权 12 形式 $\Delta$ 对应 2 维 $13$-adic Galois 表示 $\rho:\mathrm{Gal}(\bar{\mathbb Q}/\mathbb Q)\to GL_2(\mathbb Z_{13})$，满足对素数 $p\ne13$，$\mathrm{tr}\,\rho(\mathrm{Frob}_p)=\tau(p)$、$\det\rho(\mathrm{Frob}_p)=p^{11}$。于是 Hecke 关系 ($\star$) 恰是"特征多项式 $\det(1-\rho(\mathrm{Frob}_p)X)=1-\tau(p)X+p^{11}X^2$"的系数关系。猜想 1 的 $t(n)$ 满足的关系，是这族特征多项式关系在"$13^a$-精度"下的表现。

### 7.4 现代发展：Newman 问题的解决

本文 1967 年只解决了 $m=13^a$ 时问题 (3)（无限次），(4)（正密度）未解决。此后半个世纪的进展：

- **Ono（2000）**：对每个与 $6$ 互素的 $m$，$p(n)\equiv0\pmod m$ 有**正密度**（用的是 $\ell$-adic Galois 表示与 Chebotarev 密度定理）。
- **Ahlgren–Ono（2001）** 及后续工作：证明了对一大类 $m,a$，$p(n)\equiv a\pmod m$ 有正密度；最终把 Newman 的问题 (3)、(4) 在很大范围内彻底解决。
- 这些现代证明的核心，正是本文所预演的：**模方程 / Hecke 算子 / U 算子** 把一个"递推同余"提升为"分布密度"结论，只是把本文的初等工具替换成了 **Galois 表示 + 密度定理** 这一整套机器。

### 7.5 本文在方法论上的地位

1. **初等而严格**：除 Lemma 5 外，全文只用模方程、Newton 恒等式、赋值不等式，是"模形式方法初等化"的典范。
2. **"递推同余 ⇒ 无穷多个零点 ⇒ 填满陪集"** 的三段式，成为后来处理模 $\ell$ 同余的标准套路。
3. **猜想 1、2 的提出**：明确把 Hecke 乘法性"模 $13^a$ 提升"作为一个独立问题提出，预示了 $p$-adic 模形式的框架。
4. **计算与证明的互动**：作者反复强调"巧合零点"（$c(5299\cdot13^3)$、$P(13^4\cdot22655)$）靠计算获得，但证明仍严格可验证（附录 B 讨论了为何某个"看似安全"的猜想在 $N=24\cdot1530+23$ 处才首次失败）。

---

## 8. 附注与数值

**关于 $\pi$ 的最优性**：Lemma 1 的界在多数情形取等号（附录 C 表中标注 * 的少数情形例外），其"先验证明"至今（作者写论文时）不完全，是论文留下的开放问题。

**附录 C 的模方程系数表**（节选，$\beta_{1a}$）：
$$
\beta_{11}=-11\cdot13,\ \ \beta_{12}=-36\cdot13^2,\ \ \beta_{13}=-38\cdot13^3,\ \ \beta_{14}=-20\cdot13^4,\ \ \beta_{15}=-6\cdot13^5,\ \ \beta_{16}=\beta_{17}=-13^6.
$$
由此 $U\varphi=11\psi+36\cdot13\,\psi^2+38\cdot13^2\psi^3+\cdots$，并可用它验证 $L_1=\sum p(13n-7)x^n$（例如 $[x^2]$ 一侧 $=11\cdot2+468\cdot1=490=p(19)$）。

**几个可直接检验的数值事实**：
- $c(13)\equiv-\tau(1)=-1\pmod{13}$（故 Theorem 4 的陪集代表 $-1$）。
- $p(13\cdot1-7)=p(6)=11$，$p(13\cdot2-7)=p(19)=490$，$p(13\cdot3-7)=p(32)=8349$。
- $P(13^2\cdot479\,n^2)=p\big(\frac{80951n^2+1}{24}\big)=p\big(3373n^2-\frac{n^2-1}{24}\big)$。
- $(59^3\cdot13+1)/24=111247$，故定理 6 给出 (14)。

**参考文献**（论文原文所列，此处择要）：
1. J. N. O'Brien, *Some properties of partitions ...*, Ph.D. Thesis, Durham, 1966.
2. G. N. Watson, *Ramanujans Vermutung über Zerfällungsanzahlen*, J. Reine Angew. Math. **179** (1938), 97–128.
3. J. Lehner, *Divisibility properties of the Fourier coefficients of the modular invariant j(τ)*, Amer. J. Math. **71** (1949), 136–148.
4. A. O. L. Atkin, *Proof of a conjecture of Ramanujan*, Glasgow Math. J.
5. M. Newman, *Periodicity modulo m and divisibility properties of the partition function*, Trans. AMS **97** (1960), 225–236.
6. M. Newman, *Congruences for the coefficients of modular forms ...*, Canad. J. Math. **9** (1957), 549–552; **10** (1958), 577–586.
7. S. Ramanujan, *On certain arithmetical functions*, Trans. Cambridge Philos. Soc. **22** (1916), 159–184.
8. L. J. Mordell, *On Mr. Ramanujan's empirical expansions of modular functions*, Proc. Cambridge Philos. Soc. **19** (1919), 117–124.

---

> **结语**：Atkin–O'Brien 这篇 1967 年的论文是一座桥梁——它把 Ramanujan 时代关于 $p(5n+4)$ 等"算术巧合"的研究，通过模方程与 U 算子，系统化为一套"递推同余"机器，并大胆提出了 Hecke 型乘法的模 $13^a$ 提升猜想。读懂它，你就在通往 Serre、Swinnerton-Dyer 的模 $p$ 模形式理论，以及 Ono、Ahlgren 对 Newman 问题最终解答的道路上，迈出了关键一步。
