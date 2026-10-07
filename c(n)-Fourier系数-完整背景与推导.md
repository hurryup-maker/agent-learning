# $c(n)$ 的完整背景与推导：从模形式到 $j$-不变量

**—— Fourier 系数 $c(n)$ 的严谨入门讲义**

> 本讲义为《p(n) 与 c(n) 模 13 的幂的同余性质》(Atkin–O'Brien, 1967) 所需的全部背景知识。目标：让你**从零**理解 $c(n)$ 是什么、它如何被定义、如何被显式计算、以及它背后隐藏的 Hecke 结构与算术结构。所有关键公式都给出**完整推导**。

---

## 目录

1. [问题：$c(n)$ 是什么，为什么研究它](#1-问题)
2. [模群与模形式：严格基础](#2-模群与模形式)
3. [Eisenstein 级数与除数函数：完整 Fourier 展开推导](#3-eisenstein-级数)
4. [模判别式 $\Delta$ 与 Ramanujan $\tau$ 函数](#4-模判别式)
5. [$j$-不变量与 $c(n)$ 的显式计算](#5-j-不变量)
6. [Hecke 算子与 $\tau(n)$ 的乘法性：完整推导](#6-hecke-算子)
7. [$c(n)$ 与 $\tau(n)$ 的桥梁：本文要用的关键一步](#7-cn-与-taun)
8. [宏观归纳：三层结构、字典与高观点](#8-宏观归纳)
9. [附录：数值表与公式速查](#9-附录)

---

## 1. 问题

### 1.1 $c(n)$ 的定义

设 $x=e^{2\pi i\tau}$（$\tau$ 在上半平面 $\mathbb H=\{\tau:\mathrm{Im}\,\tau>0\}$）。定义
$$
j(\tau)=\frac{\big(1+240\sum_{n\ge1}\sigma_3(n)x^n\big)^3}{x\prod_{m\ge1}(1-x^m)^{24}}-744,\qquad \sigma_3(n)=\sum_{d\mid n}d^3.
$$
把 $j(\tau)$ 展开成 $x$ 的（Laurent）级数：
$$
j(\tau)=\sum_{n=-1}^{\infty}c(n)x^n.
$$
这组整数 $c(n)$ 就是本文的研究对象。前几项是
$$
c(-1)=1,\quad c(0)=0,\quad c(1)=196884,\quad c(2)=21493760,\quad c(3)=864299970,\dots
$$

> **记号约定**：通常的 Klein 模不变量 $j$ 的展开是 $j=q^{-1}+744+196884q+\cdots$（常数项为 $744$）。本文作者把 $j$ 平移了 $-744$，使其常数项为 $0$。这只影响 $c(0)$（$0$ 或 $744$），对本文关心的 $c(n),\ n\ge1$ 无影响。下文中我会在两种约定之间明确切换。

### 1.2 为什么要研究 $c(n)$？

$c(n)$ 不是"随便"定义的一串数，它同时编码了三类深刻的数学：

1. **除数函数的算术**：分子里的 $\sigma_3(n)=\sum_{d|n}d^3$；
2. **Ramanujan $\tau$ 函数的算术**：分母里的 $x\prod(1-x^m)^{24}=\sum\tau(n)x^n$ 是著名的模判别式 $\Delta$；
3. **模形式（模函数）的整体对称性**：$j$ 在 $SL_2(\mathbb Z)$ 的作用下不变。

论文的整个方法，就是利用这三种结构的相互作用，证明 $c(n)$ 模 $13$ 的幂的精细同余性质。要理解这些同余，必须先彻底理解 $c(n)$ 的"出生"过程。下面从头建立。

---

## 2. 模群与模形式

### 2.1 模群 $SL_2(\mathbb Z)$

**定义**。模群
$$
\Gamma=SL_2(\mathbb Z)=\left\{\begin{pmatrix}a&b\\c&d\end{pmatrix}:a,b,c,d\in\mathbb Z,\ ad-bc=1\right\}.
$$
它通过**分式线性变换**作用在上半平面 $\mathbb H$ 上：
$$
\gamma\cdot\tau=\frac{a\tau+b}{c\tau+d},\qquad \gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}.
$$

**生成元**。$\Gamma$ 由两个矩阵生成：
$$
S=\begin{pmatrix}0&-1\\1&0\end{pmatrix}:\ \tau\mapsto-\frac1\tau,
\qquad
T=\begin{pmatrix}1&1\\0&1\end{pmatrix}:\ \tau\mapsto\tau+1.
$$
（$S$ 是"反演"，$T$ 是"平移"。）

**基本域**。$\Gamma$ 在 $\mathbb H$ 上的基本域可取为
$$
\mathcal F=\left\{\tau\in\mathbb H:\ |\tau|\ge1,\ -\frac12\le\mathrm{Re}\,\tau\le\frac12\right\}.
$$
（边界需适当粘合。）记两个特殊点：
$$
i=(\text{阶 2 的点}),\qquad \rho=e^{2\pi i/3}=\frac{-1+i\sqrt3}{2}\ (\text{阶 3 的点}).
$$

### 2.2 模形式的定义

**定义**。设 $k\in\mathbb Z$。一个函数 $f:\mathbb H\to\mathbb C$ 称为**权 $k$ 的模形式**，若：

1. $f$ 在 $\mathbb H$ 上全纯；
2. 对一切 $\gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in\Gamma$：
$$
f\!\Big(\frac{a\tau+b}{c\tau+d}\Big)=(c\tau+d)^k f(\tau);
$$
3. $f$ 在 $\infty$ 处全纯（见下）。

**关于"在 $\infty$ 处全纯"**。由条件 2 取 $\gamma=T$，得 $f(\tau+1)=f(\tau)$，故 $f$ 有 **Fourier 展开**（$x=e^{2\pi i\tau}$）：
$$
f(\tau)=\sum_{n\in\mathbb Z}a(n)x^n.
$$
条件 3 要求 $a(n)=0$ 对一切 $n<0$，即
$$
f(\tau)=\sum_{n=0}^{\infty}a(n)x^n.
$$
此时称 $f$ 在 $\infty$ 全纯，并记 $f(\infty)=a(0)$。

若还满足 $a(0)=0$，则称 $f$ 为**尖点形式**（cusp form）。

**权 $k$ 模形式空间**记作 $M_k$，尖点形式空间记作 $S_k$。二者都是有限维复向量空间。

### 2.3 维数公式（陈述）

对偶数 $k\ge4$：
$$
\dim M_k=\begin{cases}\lfloor k/12\rfloor,& k\equiv2\pmod{12},\\[2pt]
\lfloor k/12\rfloor+1,& k\not\equiv2\pmod{12},\end{cases}
\qquad
\dim S_k=\dim M_k-1\ \ (k\ge4).
$$
我们需要的两个特例：
$$
\dim M_4=1,\quad \dim M_6=1,\quad \dim M_{12}=2,\quad \dim S_{12}=1.
$$
**$S_{12}$ 是一维的**——这将是第 4、6 节所有论证的支点。

### 2.4 零点公式（valence formula）

**定理（零点公式）**。设 $f\ne0$ 是权 $k$ 的模形式。则
$$
\operatorname{ord}_\infty(f)+\frac12\operatorname{ord}_i(f)+\frac13\operatorname{ord}_\rho(f)+\sum_{\substack{P\in\mathcal F\\P\ne i,\rho}}\operatorname{ord}_P(f)=\frac{k}{12}.
\tag{VF}
$$
这里 $\operatorname{ord}_P(f)$ 是 $f$ 在 $P$ 的零点阶数；$i,\rho$ 的系数 $1/2,1/3$ 来自它们的稳定子阶数 $2,3$。

**证明思路（围道积分）**。考虑 $f'/f$ 沿基本域边界 $\partial\mathcal F$ 的积分 $\frac{1}{2\pi i}\oint\frac{f'}f\,d\tau$。由留数定理，它等于 $\mathcal F$ 内零点阶数（加权）之和。另一方面，边界各段通过 $S,T$ 的变换两两抵消（利用 $f(\gamma\tau)=(c\tau+d)^kf(\tau)$ 取对数导数产生的项 $(c\tau+d)'/(c\tau+d)=c/(c\tau+d)$），最终算出弧段贡献 $-k/12$。整理即得 (VF)。$\blacksquare$

**关键应用**：任何非零尖点形式 $f\in S_k$ 在 $\infty$ 处 $\operatorname{ord}_\infty(f)\ge1$。由 (VF)，$k/12\ge1$ 且等号当且仅当 $f$ 在 $\mathbb H$ 内无零点。特别地，权 $12$ 的尖点形式 $\Delta$（见下）满足 $12/12=1$，故 **$\Delta$ 在 $\mathbb H$ 上无零点**。这是 $j=\text{（某物）}/\Delta$ 有意义的根本原因。

---

## 3. Eisenstein 级数

### 3.1 定义

对偶数 $k\ge4$，定义 **Eisenstein 级数**
$$
G_k(\tau)=\sum_{\substack{(m,n)\in\mathbb Z^2\\(m,n)\ne(0,0)}}\frac1{(m\tau+n)^k}.
$$
（级数绝对收敛，因 $k\ge4$。）规范化为
$$
E_k(\tau)=\frac{G_k(\tau)}{2\zeta(k)},\qquad \zeta(k)=\sum_{n\ge1}\frac1{n^k}.
$$

### 3.2 两个引理（$E_k$ 是模形式 + Fourier 展开）

**引理 A（模性）**。$G_k$（从而 $E_k$）是权 $k$ 的模形式。

**证明**。对 $\gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}$，
$$
G_k(\gamma\tau)=\sum_{(m,n)\ne0}\Big(\frac{am+b}{c\tau+d}...\Big)^{-k}
$$
具体地，$G_k(\gamma\tau)=\sum_{(m,n)\ne0}\big(\frac{a\cdot\frac{a\tau+b}{c\tau+d}+b}{...}\big)^{-k}$。把 $(m,n)\mapsto(m',n')=(am+cn,\ bm+dn)$ 是 $\mathbb Z^2$ 上的双射（因 $\det\gamma=1$），且
$$
m\cdot\frac{a\tau+b}{c\tau+d}+n=\frac{m(a\tau+b)+n(c\tau+d)}{c\tau+d}=\frac{(am+cn)\tau+(bm+dn)}{c\tau+d}=\frac{m'\tau+n'}{c\tau+d}.
$$
故
$$
G_k(\gamma\tau)=(c\tau+d)^k\sum_{(m',n')\ne0}\frac1{(m'\tau+n')^k}=(c\tau+d)^kG_k(\tau).
$$
在 $\infty$ 处全纯由下面引理 B 的展开看出。$\blacksquare$

### 3.3 Lipschitz 公式（完整推导）

**引理 B（Lipschitz）**。对 $k\ge2$、$\mathrm{Im}\,z>0$：
$$
\sum_{n\in\mathbb Z}\frac1{(z+n)^k}=\frac{(-1)^k(2\pi i)^k}{(k-1)!}\sum_{r\ge1}r^{k-1}e^{2\pi irz}.
\tag{Lip}
$$

**证明**。分两步。

**第 1 步：余切的部分分式展开。**
$$
\pi\cot(\pi z)=\frac1z+\sum_{n\ge1}\Big(\frac1{z-n}+\frac1{z+n}\Big)=\sum_{n\in\mathbb Z}\frac1{z+n}
\tag{1}
$$
（对称求和，避免发散）。这是经典的 $\cot$ 部分分式公式（可由 $\sin$ 的 Weierstrass 无穷乘积 $\sin(\pi z)=\pi z\prod_{n\ge1}(1-z^2/n^2)$ 取对数导数得到）。

**第 2 步：把 $\cot$ 写成 $q$ 级数。**
$$
\cot(\pi z)=\frac{\cos(\pi z)}{\sin(\pi z)}
=\frac{i(e^{i\pi z}+e^{-i\pi z})}{e^{i\pi z}-e^{-i\pi z}}
=i\,\frac{e^{2\pi iz}+1}{e^{2\pi iz}-1}.
$$
令 $q=e^{2\pi iz}$（$\mathrm{Im}\,z>0$ 保证 $|q|<1$），则
$$
\cot(\pi z)=i\,\frac{q+1}{q-1}=-i\,\frac{1+q}{1-q}
=-i\,(1+q)(1+q+q^2+\cdots)=-i-2i\sum_{r\ge1}q^r.
$$
于是
$$
\pi\cot(\pi z)=-\pi i-2\pi i\sum_{r\ge1}e^{2\pi irz}.
\tag{2}
$$

**第 3 步：对 (1) 微分 $k-1$ 次。**
对 (1) 两边取 $\dfrac{d^{k-1}}{dz^{k-1}}$：
$$
\text{左端}=\sum_{n\in\mathbb Z}\frac{(-1)(-2)\cdots(-(k-1))}{(z+n)^k}
=\sum_{n\in\mathbb Z}\frac{(-1)^{k-1}(k-1)!}{(z+n)^k}.
$$
对 (2) 右端取 $\dfrac{d^{k-1}}{dz^{k-1}}$（常数项 $-\pi i$ 消失）：
$$
\frac{d^{k-1}}{dz^{k-1}}\big[{-}2\pi i e^{2\pi irz}\big]
={-}2\pi i\,(2\pi ir)^{k-1}e^{2\pi irz}.
$$
故
$$
\sum_{n\in\mathbb Z}\frac{(-1)^{k-1}(k-1)!}{(z+n)^k}
={-}2\pi i\sum_{r\ge1}(2\pi ir)^{k-1}e^{2\pi irz}
={-}(2\pi i)^k\sum_{r\ge1}r^{k-1}e^{2\pi irz}.
$$
两边除以 $(-1)^{k-1}(k-1)!$，注意 $(-1)^{k-1}\cdot(-1)^{k-1}=1$：
$$
\sum_{n\in\mathbb Z}\frac1{(z+n)^k}
=\frac{(-1)^{k-1}(k-1)!}{(-1)^{k-1}(k-1)!}\cdot\frac{-(2\pi i)^k}{(-1)^{k-1}(k-1)!}\sum r^{k-1}e^{2\pi irz}
=\frac{(-1)^k(2\pi i)^k}{(k-1)!}\sum_{r\ge1}r^{k-1}e^{2\pi irz}.
$$
（因 $-(2\pi i)^k/(-1)^{k-1}=(-1)^k(2\pi i)^k$。）即 (Lip)。$\blacksquare$

### 3.4 $G_k$ 与 $E_k$ 的 Fourier 展开（完整推导）

把 $G_k$ 的和拆成 $m=0$ 与 $m\ne0$：
$$
G_k(\tau)=\underbrace{\sum_{n\ne0}\frac1{n^k}}_{m=0\ \text{项}}+\sum_{m\ne0}\sum_{n\in\mathbb Z}\frac1{(m\tau+n)^k}.
$$
第一项 $=\sum_{n\ne0}n^{-k}=2\zeta(k)$（$k$ 偶）。第二项对 $m\ne0$ 用 Lipschitz（$z=m\tau$）：
$$
\sum_{n\in\mathbb Z}\frac1{(m\tau+n)^k}=\frac{(-1)^k(2\pi i)^k}{(k-1)!}\sum_{r\ge1}r^{k-1}e^{2\pi irm\tau}.
$$
现在 $k=2\ell$（偶）。注意 $(2\pi i)^{2\ell}=(2\pi)^{2\ell}(-1)^\ell$，$(-1)^{2\ell}=1$。故
$$
\sum_{n\in\mathbb Z}\frac1{(m\tau+n)^{2\ell}}=\frac{(2\pi)^{2\ell}(-1)^\ell}{(2\ell-1)!}\sum_{r\ge1}r^{2\ell-1}e^{2\pi irm\tau}.
$$
把 $m\ne0$ 的两边（$m>0$ 与 $m<0$ 贡献相同，因 $e^{2\pi irm\tau}$ 只依赖 $|m|$）合并：
$$
G_{2\ell}(\tau)=2\zeta(2\ell)+2\,\frac{(2\pi)^{2\ell}(-1)^\ell}{(2\ell-1)!}\sum_{m\ge1}\sum_{r\ge1}r^{2\ell-1}e^{2\pi irm\tau}.
$$
令 $n=mr$，则 $\sum_{m,r\ge1}r^{2\ell-1}e^{2\pi irm\tau}=\sum_{n\ge1}\big(\sum_{r\mid n}r^{2\ell-1}\big)q^n=\sum_{n\ge1}\sigma_{2\ell-1}(n)q^n$（$q=e^{2\pi i\tau}$）。得
$$
\boxed{G_{2\ell}(\tau)=2\zeta(2\ell)+\frac{2(-1)^\ell(2\pi)^{2\ell}}{(2\ell-1)!}\sum_{n\ge1}\sigma_{2\ell-1}(n)q^n.}
\tag{3}
$$

除以 $2\zeta(2\ell)$，并用 $\zeta(2\ell)=\dfrac{(-1)^{\ell+1}(2\pi)^{2\ell}B_{2\ell}}{2(2\ell)!}$（Euler，$B$ 为 Bernoulli 数）：
$$
E_{2\ell}(\tau)=1+\frac{(-1)^\ell(2\pi)^{2\ell}}{(2\ell-1)!\zeta(2\ell)}\sum\sigma_{2\ell-1}(n)q^n.
$$
计算系数（这是本推导唯一的代数耐心处）：
$$
\frac{(-1)^\ell(2\pi)^{2\ell}}{(2\ell-1)!\zeta(2\ell)}
=\frac{(-1)^\ell(2\pi)^{2\ell}}{(2\ell-1)!}\cdot\frac{(-1)^{\ell+1}2(2\ell)!}{(2\pi)^{2\ell}B_{2\ell}}
=\frac{2(2\ell)!}{(2\ell-1)!B_{2\ell}}\cdot(-1)^{2\ell+1}
=-\frac{2(2\ell)!}{(2\ell-1)!B_{2\ell}}
=-\frac{4\ell}{B_{2\ell}}.
$$
故得到**核心公式**：
$$
\boxed{E_{2\ell}(\tau)=1-\frac{4\ell}{B_{2\ell}}\sum_{n\ge1}\sigma_{2\ell-1}(n)q^n.}
\tag{4}
$$

### 3.5 $E_4$ 与 $E_6$ 的显式

代入 Bernoulli 数 $B_4=-\frac1{30}$、$B_6=\frac1{42}$：

$$
E_4(\tau)=1-\frac{8}{B_4}\sum\sigma_3(n)q^n=1+240\sum_{n\ge1}\sigma_3(n)q^n,
\tag{5}
$$

$$
E_6(\tau)=1-\frac{12}{B_6}\sum\sigma_5(n)q^n=1-504\sum_{n\ge1}\sigma_5(n)q^n.
\tag{6}
$$

展开前几项：
$$
E_4=1+240q+2160q^2+6720q^3+17520q^4+\cdots,
$$
$$
E_6=1-504q-16632q^2-122976q^3-\cdots.
$$

> **要点**：这就是 $c(n)$ 定义式分子里那个 $1+240\sum\sigma_3(n)q^n$ 的来源——它是 $E_4$。从格点求和一路推到除数函数 $\sigma_3(n)$，是理解 $c(n)$ 算术内涵的第一块拼图。

---

## 4. 模判别式 $\Delta$

### 4.1 定义

定义
$$
\Delta(\tau)=\frac{E_4^3-E_6^2}{1728}.
\tag{7}
$$
为什么是 $1728=12^3$？因为 $E_4\equiv E_6\equiv1\pmod q$，故
$$
E_4^3-E_6^2=(1+720q+\cdots)-(1-1008q+\cdots)=1728q+\cdots,
$$
除以 $1728$ 使 $\Delta$ 的 $q$-系数为整数（且首项为 $1$）。由 (7)，$\Delta$ 是权 $12$ 的模形式，且 $\Delta(\infty)=0$，故 $\Delta\in S_{12}$。

### 4.2 $\Delta=\eta^{24}$：完整证明

**定理**。$\Delta(\tau)=x\prod_{m\ge1}(1-x^m)^{24}=x\,f^{24}(x)$，其中 $x=e^{2\pi i\tau}$。

**证明**。设 $g(\tau)=x\prod(1-x^m)^{24}$。因 $\eta(\tau)=x^{1/24}\prod(1-x^m)$ 是权 $\tfrac12$ 的模形式（带乘子系），故 $\eta^{24}=x\prod(1-x^m)^{24}=g$ 是权 $12$ 的模形式；且 $g=x(1-x)^{24}\cdots=x-24x^2+\cdots$ 无常数项，是尖点形式，即 $g\in S_{12}$。

由维数公式，$\dim S_{12}=1$。故 $\Delta$ 与 $g$ 线性相关：存在常数 $\lambda$ 使 $\Delta=\lambda g$。比较 $q$ 的系数：$\Delta$ 的首项为 $1728q/1728=q$（因 $E_4^3-E_6^2=1728q+\cdots$），$g$ 的首项为 $q$。故 $\lambda=1$，即 $\Delta=g=x\prod(1-x^m)^{24}$。$\blacksquare$

### 4.3 Ramanujan $\tau$ 函数

展开
$$
\Delta(\tau)=x\prod_{m\ge1}(1-x^m)^{24}=\sum_{n\ge1}\tau(n)x^n.
$$
前几项（可由 Euler 五边形数定理逐项算得）：
$$
\begin{array}{c|rrrrrrrr}
n&1&2&3&4&5&6&7&8\\\hline
\tau(n)&1&-24&252&-1472&4830&-6048&-16744&84480
\end{array}
$$
$\tau(n)$ 就是 **Ramanujan $\tau$ 函数**。Ramanujan 在 1916 年对它提出三个著名猜想：乘法性（见 §6）、$\tau(n)\le d(n)n^{11/2}$ 的增长估计（"Ramanujan 猜想"，Deligne 1974 证明）、以及 $\tau(n)\equiv\sigma_{11}(n)\pmod{691}$ 等模同余。

> 本文中 $\Delta=x f^{24}(x)$ 恰好是论文里 $\sum r(n)x^n=x f^{24}(x)$ 的定义，即 $r(n)=\tau(n)$。

---

## 5. $j$-不变量

### 5.1 定义与良好定义性

**定义**。
$$
j(\tau)=\frac{E_4^3(\tau)}{\Delta(\tau)}.
\tag{8}
$$

**定理（$j$ 是权 $0$ 的模函数）**。$j$ 在 $\mathbb H$ 上全纯，且对一切 $\gamma\in\Gamma$，$j(\gamma\tau)=j(\tau)$。

**证明**。全纯性：$E_4^3$ 全纯，而由零点公式，$\Delta$ 在 $\mathbb H$ 上无零点（§2.4），故商全纯。不变性：$E_4$ 权 $4$、$\Delta$ 权 $12$，故
$$
j(\gamma\tau)=\frac{E_4^3(\gamma\tau)}{\Delta(\gamma\tau)}
=\frac{\big((c\tau+d)^4E_4(\tau)\big)^3}{(c\tau+d)^{12}\Delta(\tau)}
=\frac{(c\tau+d)^{12}E_4^3}{(c\tau+d)^{12}\Delta}=j(\tau).
\blacksquare
$$

（"权 $0$ 模函数"就是"$\Gamma$-不变的全纯函数，允许在 $\infty$ 有极点"。）

### 5.2 为什么是 $j$？——它生成整个模函数域

**定理**。$SL_2(\mathbb Z)$ 的模函数域恰为 $\mathbb C(j)$：任何在 $\mathbb H$ 全纯、$\Gamma$-不变、且在 $\infty$ 至多有极点的亚纯函数，都是 $j$ 的有理函数。

（证明思路：若 $f$ 是权 $0$ 模函数，选常数 $A$ 使 $f-A$ 的零点抵消 $j$ 在 $i,\rho$ 等的零点，再用零点公式说明 $j$ 的极点阶数恰为 $1$，从而 $j$ 是"分母 $\Delta$、分子 $E_4^3$"的最小组合。$j$ 在 $\infty$ 有**单极点**。）

### 5.3 $c(n)$ 的显式计算（完整推导前几项）

$j=E_4^3/\Delta$。用 $\Delta=x f^{24}(x)$ 与 $E_4=1+240\sum\sigma_3(n)x^n$，这正是论文里的定义式（去掉 $-744$ 平移）：
$$
j=\frac{E_4^3}{\Delta}=\frac{(1+240\sum\sigma_3(n)x^n)^3}{x f^{24}(x)}.
$$

**算 $1/\Delta$**。因 $f^{24}(x)=\prod(1-x^m)^{24}$，其倒数
$$
\frac1{\Delta}=\frac1{x\,f^{24}(x)}=\frac1x\prod_{m\ge1}\frac1{(1-x^m)^{24}}
=\frac1x\sum_{n\ge0}p_{24}(n)x^n,
$$
其中 $p_{24}(n)$ 是"$24$ 色分拆数"（把 $n$ 分拆成若干部分、每部分染上 $24$ 种颜色之一的方法数）。于是
$$
\frac1{\Delta}=x^{-1}+24+324x+3200x^2+25650x^3+\cdots
$$
（由 $p_{24}(0)=1,\ p_{24}(1)=24,\ p_{24}(2)=324,\ p_{24}(3)=3200,\dots$）。

**乘上 $E_4^3$**。$E_4=1+240x+2160x^2+6720x^3+\cdots$，故
$$
E_4^3=(1+240x+2160x^2+\cdots)^3=1+720x+179280x^2+14800800x^3+\cdots.
$$
（验证：$q$ 系数 $3\cdot240=720$；$q^2$ 系数 $3\cdot2160+3\cdot240^2=6480+172800=179280$。）

**相乘**：
$$
j=E_4^3\cdot\frac1{\Delta}=(1+720x+179280x^2+\cdots)(x^{-1}+24+324x+3200x^2+\cdots).
$$
逐项收集：
- $x^{-1}$：$1$；
- $x^0$：$24+720=744$；
- $x^1$：$324+720\cdot24+179280=324+17280+179280=196884$；
- $x^2$：$3200+720\cdot324+179280\cdot24+14800800=21493760$。

故
$$
\boxed{j=x^{-1}+744+196884x+21493760x^2+864299970x^3+\cdots}
\tag{9}
$$
与标准值完全一致。所以
$$
c(-1)=1,\ c(0)=744,\ c(1)=196884,\ c(2)=21493760,\ c(3)=864299970,\dots
$$

### 5.4 本文的归一化（$-744$）

论文写 $j=(1+240\sum\sigma_3x^n)^3/(x f^{24}(x))-744$，即
$$
j_{\text{论文}}=j-744=x^{-1}+0+196884x+21493760x^2+\cdots.
$$
所以在论文的记号下 $c(0)=0$。**本文所有同余都只涉及 $n\ge1$，两种归一化无差别**。

---

## 6. Hecke 算子

### 6.1 定义

对权 $k$ 模形式 $f=\sum_{n\ge0}a(n)q^n$ 与正整数 $m$，定义 **Hecke 算子 $T_m$**：
$$
T_m f(\tau)=m^{k-1}\sum_{\substack{ad=m\\ b=0,\dots,d-1}}\frac1{d^k}f\!\Big(\frac{a\tau+b}{d}\Big).
\tag{10}
$$
（这是原始定义；等价地，$T_m$ 由 $f$ 的 $q$ 系数公式给出，见下。）

**素数情形 $T_p$ 的 $q$ 展开作用**（我们只需要这个）：设 $f=\sum a(n)q^n$，则
$$
T_p f(\tau)=\sum_{n\ge0}\big(a(pn)+p^{k-1}a(n/p)\big)q^n,
\tag{11}
$$
其中约定 $n/p$ 非整数时 $a(n/p)=0$。**关键**：$T_p f$ 仍是权 $k$ 模形式。

**与本文 $U$ 算子的关系**。本文的 $U$（取 $13$ 的倍数项）
$$
Uf=\sum a(pn)q^n
$$
正是 (11) 的第一项。因此 $T_p=U+p^{k-1}V_p$，其中 $V_p f=\sum a(n/p)q^n=f(q^p)$。这说明论文里的 $U$ 就是 Hecke 算子在 $p\mid$ 时的"主导项"。

### 6.2 $\Delta$ 是特征形式（完整推导）

**定理**。$T_m\Delta=\tau(m)\Delta$ 对一切 $m\ge1$。

**证明**。因 $\dim S_{12}=1$ 且 $T_m$ 保持 $S_{12}$，故存在 $\lambda_m$ 使 $T_m\Delta=\lambda_m\Delta$。比较 $q$ 的系数：由 (11)（$k=12$，$k-1=11$），$T_m\Delta$ 的 $q$ 系数为
$$
a(1\cdot m)+m^{11}a(1/m)=a(m)+0=\tau(m).
$$
而 $\lambda_m\Delta$ 的 $q$ 系数为 $\lambda_m\tau(1)=\lambda_m$。故 $\lambda_m=\tau(m)$。$\blacksquare$

### 6.3 $\tau$ 的乘法性（完整推导）

**定理（Ramanujan）**。$\tau$ 是**乘性**的：$\tau(mn)=\tau(m)\tau(n)$ 当 $(m,n)=1$；且对素数 $p$、任意 $m$：
$$
\tau(pm)=\tau(p)\tau(m)-p^{11}\tau(m/p).
\tag{12}
$$

**证明**。由 $T_p\Delta=\tau(p)\Delta$ 与 (11)，比较两边 $q^m$ 系数：
$$
\tau(pm)+p^{11}\tau(m/p)=\tau(p)\tau(m),
$$
移项即得 (12)。（互素情形由 $T_mT_n=T_{mn}$（$(m,n)=1$）与特征形式性质 $\lambda_{mn}=\lambda_m\lambda_n$ 推出。）$\blacksquare$

**例**：$\tau(6)=\tau(2\cdot3)=\tau(2)\tau(3)=(-24)(252)=-6048$（与上表一致）；$\tau(4)=\tau(2)\tau(2)-2^{11}\tau(1)=576-2048=-1472$。$\blacksquare$

> **这正是本文猜想 1 的"原型"**。论文把 (12) 中 $\tau$ 换成 $t(n)\equiv c(13^a n)/c(13^a)\pmod{13^a}$，猜想同样的乘法关系模 $13^a$ 成立——这就是"把 Hecke 特征值的算术提升到 $13^a$ 精度"。

---

## 7. $c(n)$ 与 $\tau(n)$ 的桥梁

### 7.1 模 $13$ 的结构

本文最关键的一步是建立 $c(13n)$ 与 $\tau(n)$ 的关系。其根源如下。

由 §5，$j-744=E_4^3/\Delta-744$。模 $13$ 时，$E_4\equiv1+240\sum\sigma_3(n)q^n\pmod{13}$。而 $\sigma_3(n)\equiv\sigma_1(n)=\sum_{d|n}d\pmod{13}$（因 $d^3\equiv d$ 当 $(d,13)=1$，且 $13$ 的幂贡献 $0$），$\sigma_1(n)\equiv0\pmod{13}$ 当 $13\mid n$。经过系统化整理（这正是论文 §2–§3 的 Lemma 4、5 所做），得到：

$$
\sum_{n\ge1}c(13n)x^n=-\psi(x)+13^2U^2\psi(x),\qquad \psi(x)=x\frac{f^2(x^{13})}{f^2(x)}.
\tag{L5}
$$

模 $13$ 取首项，再用 **Frobenius 同态**（"新手梦"）$f(x^{13})\equiv f(x)^{13}\pmod{13}$：
$$
\psi(x)=x\frac{f^2(x^{13})}{f^2(x)}\equiv x\frac{f^{26}(x)}{f^2(x)}=x f^{24}(x)=\Delta(x)\pmod{13}.
$$
于是
$$
\boxed{c(13n)\equiv-\tau(n)\pmod{13}.}
\tag{13}
$$
这就是 (L5) 的第一个推论，也是 Newman 证明 $c(13^2n)\equiv8c(13n)\pmod{13}$ 的出发点。

### 7.2 一句话总结桥梁

$c(n)$ 的"第 $13$ 层"（$c(13n)$）在模 $13$ 下**恰好等于**（差一个符号）Ramanujan $\tau$ 函数。而 $\tau$ 是 Hecke 特征形式（§6），有漂亮的乘法结构。于是 $\tau$ 的 Hecke 算术"渗透"到 $c(n)$ 上——论文的猜想 1 就是要把这种渗透提升到 $13$ 的任意幂。

---

## 8. 宏观归纳

### 8.1 三层结构

把 $c(n)$ 的知识组织成三层，由浅入深：

| 层次 | 对象 | 结构 |
|---|---|---|
| **0. 形式层** | $E_4=1+240\sum\sigma_3(n)q^n$，$\Delta=q\prod(1-q^m)^{24}$，$j=E_4^3/\Delta$ | 纯幂级数运算，可逐项算 $c(n)$ |
| **1. 模层** | $E_4\in M_4$，$\Delta\in S_{12}$，$j$ 是模函数 | 模性 + 维数/零点公式 ⇒ $\Delta=\eta^{24}$、$j$ 生成模函数域 |
| **2. 算术层** | $\tau(n)$ 是 Hecke 特征值，$\sigma_3(n)$ 是除数函数 | $T_p\Delta=\tau(p)\Delta$ ⇒ $\tau$ 乘性；模 $13$：$c(13n)\equiv-\tau(n)$ |

论文的工作在第 2 层：**用 Hecke 算术 + 模方程，把 $c(13^a n)$ 的模 $13^a$ 行为完全刻划**。

### 8.2 一张"字典"

$$
\begin{array}{c|ccc}
& \text{权} & \text{原型} & \text{系数}\\
\hline
\text{Eisenstein }E_4 & 4 & \text{格点求和} & \sigma_3(n)\ (\text{除数函数})\\
\text{判别式 }\Delta & 12\ (\text{尖点}) & \eta^{24} & \tau(n)\ (\text{Hecke 特征值})\\
j\text{-不变量} & 0\ (\text{极点在}\infty) & E_4^3/\Delta & c(n)\ (\text{二者之商})
\end{array}
$$

$c(n)$ 是"除数函数算术"与"$\tau$ 算术"的**商**，这正是它既有复杂的分布、又隐藏着可刻划的同余结构的原因。

### 8.3 高观点速览（不展开，仅定位）

1. **Hecke 特征值视角**：$c(13n)\equiv-\tau(n)\pmod{13}$ 意味着 $c(n)$ 的"第 13 层"继承了权 12 特征形式 $\Delta$ 的全部 Hecke 算术。
2. **Galois 表示视角**（Deligne 1974）：$\tau(p)$ 是某个 2 维 $\ell$-adic 表示的 Frobenius 迹，$\det=p^{11}$。于是 (12) 是特征多项式 $\det(1-\rho(\mathrm{Frob}_p)X)=1-\tau(p)X+p^{11}X^2$ 的系数关系。
3. **Moonshine（趣闻）**：$c(1)=196884=1+196883$，其中 $196883$ 是"魔群"（Monster）最小非平凡不可约表示的维数。这是 McKay 观察、后被 Borcherds 证明的"怪兽月光"的起点。
4. **增长性**：$c(n)\sim\dfrac{e^{4\pi\sqrt n}}{\sqrt2\,n^{3/4}}$（由 $j$ 在 $\infty$ 的主项 $q^{-1}$ 与 Hardy–Ramanujan 型方法得到），故 $c(n)$ 增长极快，但其模 $13$ 行为却高度规律——这正是模形式"局部（模 $\ell$）规律、全局（实值）爆炸"的典型现象。

---

## 9. 附录

### 9.1 常用公式速查

- **Euler 五边形数定理**：$f(x)=\prod_{m\ge1}(1-x^m)=\sum_{k=-\infty}^{\infty}(-1)^kx^{k(3k-1)/2}$。
- **Frobenius 同态**：对素数 $p$，$(1-x)^{p}\equiv1-x^p\pmod p$；故 $f(x^{p})=\prod(1-x^{pm})\equiv\prod(1-x^m)^{p}=f(x)^p\pmod p$。
- **除数函数**：$\sigma_k(n)=\sum_{d\mid n}d^k$，乘性；对 $n=\prod p_i^{e_i}$，$\sigma_k(n)=\prod\frac{p_i^{k(e_i+1)}-1}{p_i^k-1}$。
- **Bernoulli 数**：$B_2=\tfrac16,\ B_4=-\tfrac1{30},\ B_6=\tfrac1{42},\ B_8=-\tfrac1{30},\ B_{10}=\tfrac5{66},\ B_{12}=-\tfrac{691}{2730}$。
- **$\zeta$ 值**：$\zeta(2)=\tfrac{\pi^2}{6},\ \zeta(4)=\tfrac{\pi^4}{90},\ \zeta(6)=\tfrac{\pi^6}{945}$。

### 9.2 前几项数值表

$$
\begin{array}{c|cccccccc}
n & -1 & 0 & 1 & 2 & 3 & 4 & 5\\\hline
c(n)\ (\text{标准 }j) & 1 & 744 & 196884 & 21493760 & 864299970 & 20245856256 & 333202640600\\
c(n)\ (\text{论文 }j-744) & 1 & 0 & 196884 & 21493760 & 864299970 & 20245856256 & 333202640600\\
\tau(n) & & & 1 & -24 & 252 & -1472 & 4830\\
\sigma_3(n) & & & 1 & 9 & 28 & 73 & 126
\end{array}
$$

### 9.3 一个可直接验证的模 13 关系

$c(13)\equiv-\tau(1)\equiv-1\equiv12\pmod{13}$（论文用它作为定理 4 的陪集代表）。更一般地，由 (13)，$c(13n)\equiv-\tau(n)\pmod{13}$ 对一切 $n$ 成立，而 $\tau$ 的乘法性 (12) 给出了 $c(13n)$ 的"乘法结构"，这正是本文一切结论的种子。

---

> **结语**：$c(n)$ 的三个来源——$E_4$ 的除数函数、$\Delta$ 的 $\tau$ 函数、$j$ 的模不变性——共同决定了它的算术行为。本文作者抓住"$c(13n)\equiv-\tau(n)\pmod{13}$"这座桥，把 Hecke 理论搬上了 $13^a$ 的阶梯。现在你有了完整的背景，可以回到主讲义，从容地读每一个引理和定理了。
