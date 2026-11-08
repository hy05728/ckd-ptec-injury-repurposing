# 推送到 GitHub —— 两条命令

仓库已在本机初始化并提交完毕（commit `dbc35ff`，42 个文件，1.99 MB）。
你只需要：**建一个空仓库 → 关联 → 推送**。

目标地址：`https://github.com/hy05728/ckd-ptec-injury-repurposing`

---

## 第 0 步：在 GitHub 网页建空仓库（必做，顺序不能颠倒）

1. 打开 https://github.com/new
2. **Repository name** 填：`ckd-ptec-injury-repurposing`（必须一字不差，否则要改文章里的 URL）
3. 可见性选 **Public**（文章写了公开可获取）
4. **下面三个勾选框全部不要勾**：Add a README / Add .gitignore / Choose a license
   —— 勾了会产生初始提交，push 时会冲突
5. Create repository

---

## 第 1 条：关联远程

```bash
cd /e/workbuddy/2026-11-08-00-03-12/repo_github && git remote add origin https://github.com/hy05728/ckd-ptec-injury-repurposing.git
```

## 第 2 条：推送

```bash
git push -u origin main
```

---

## 推送时会要求输入账号密码 —— 密码必须用 Token，不是登录密码

GitHub 从 2021 年起不接受账号密码。弹出的密码框要填 **Personal Access Token**：

1. 打开 https://github.com/settings/tokens
2. **Generate new token** → **Generate new token (classic)**
3. Note 随便填（如 `ckd-repo-push`），Expiration 选 90 days
4. 权限只勾 **`repo`** 这一整块，其它不用管
5. Generate token，**立刻复制那串 `ghp_...`**（只显示一次）
6. 回终端：Username 填 `hy05728`，Password 粘贴 token

---

## 验证

打开 `https://github.com/hy05728/ckd-ptec-injury-repurposing`，应能看到 README 渲染、
`scripts/` 24 个脚本、`results/` 12 个结果表。

**这一步必须做成**：文章正文 Data and code availability 写了这个 URL，审稿人会点开。

---

## 想改提交邮箱？

当前用的是 `jessis0530@163.com`。若你 GitHub 账号绑的是别的邮箱：

```bash
cd /e/workbuddy/2026-11-08-00-03-12/repo_github && git config user.email "你的邮箱" && git commit --amend --reset-author --no-edit
```

---

## 常见报错

| 报错 | 原因 | 解决 |
|---|---|---|
| `remote origin already exists` | 之前加过 | `git remote set-url origin https://github.com/hy05728/ckd-ptec-injury-repurposing.git` |
| `failed to push some refs` / `rejected` | 建仓库时勾了 README，远程已有提交 | `git pull --rebase origin main` 再 `git push -u origin main` |
| `Authentication failed` | 密码填了登录密码 | 用 Personal Access Token（见上） |
| `Support for password authentication was removed` | 同上 | 同上 |
| `fatal: not a git repository` | 不在 repo_github 目录里 | `cd /e/workbuddy/2026-11-08-00-03-12/repo_github` |

---

## 仓库内容

| 路径 | 内容 |
|---|---|
| `README.md` | 数据来源表、环境、19 步流水线、关键结果、复现性说明 |
| `LICENSE` | MIT，署名 Yue Hu |
| `scripts/` | 24 个脚本：下载 → 转码 → 合并 → QC → 签名 → 筛选 → 对接 → 基准 → 富集 → 稳健性 → 出图 → 自检 |
| `md/` | AutoDL 上的 MD + MM/PBSA 完整工具链（含中文操作手册，暂未使用）|
| `results/` | 12 个关键结果表（签名、排名、对接、置换检验、稳健性）|
| `.gitignore` / `.gitattributes` | 排除大文件；统一 LF 换行，避免 Linux 上脚本失效 |

**未包含**：原始数据（1.5 GB GEO 压缩包、1.57 GB CELLxGENE h5ad），
由 `scripts/00_download.py` 下载。这是标准做法，也符合 `.gitignore` 的设定。
