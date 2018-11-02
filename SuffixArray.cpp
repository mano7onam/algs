#pragma GCC optimize ("O3")
#pragma GCC target ("sse4")
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;
#define mp make_pair

void print_vector(string name, vector<int> v) {
    cout << name << ": ";
    for (int i = 0; i < v.size(); ++i) {
        cout << v[i] << ' ';
    }
    cout << endl;
}

template<class T> struct SegTree {
    T INF = 1e9 + 7;
    vector<T> arr;
    vector<T> tree;

    void build(int v, int tl, int tr) {
        if (tl == tr) {
            tree[v] = arr[tl];
            return;
        }
        int tm = (tl + tr) / 2;
        build(2 * v, tl, tm);
        build(2 * v + 1, tm + 1, tr);
    }

    SegTree(vector<T> arr) : arr(arr) {
        tree.assign(arr.size() * 4, 0);
        build(1, 0, arr.size() - 1);
    }

    void update(int v, int tl, int tr, int pos, T val) {
        if (tl == tr) {
            tree[v] = min(tree[v], val);
            return;
        }
        int tm = (tl + tr) / 2;
        if (pos <= tm) {
            update(2 * v, tl, tm, pos, val);
        } else {
            update(2 * v + 1, tm + 1, tr, pos, val);
        }
        tree[v] = min(tree[2 * v], tree[2 * v + 1]);
    }

    T get(int v, int tl, int tr, int l, int r) {
        if (l == tl && r == tr) {
            return tree[v];
        }
        int tm = (tl + tr) / 2;
        T res = INF;
        if (l <= tm) {
            res = min(res, get(2 * v, tl, tm, l, min(r, tm)));
        }
        if (r > tm) {
            res = min(res, get(2 * v + 1, tm + 1, tr, max(l, tm + 1), r));
        }
        return res;
    }

    T get(int l, int r) {
        return get(1, 0, arr.size() - 1, l, r);
    }
};

struct SuffixArray {
    int N;
    vector<int> idx;
    string str;
    
    // 34 37 41 37 -> 1, 2, 3, 2
    void compress(vector<int>& v) {
        vector<int> V = v; 
	    sort(V.begin(), V.end()); 
	    V.erase(unique(V.begin(), V.end()), V.end());
        for (int& i: v) {
	        i = lower_bound(V.begin(), V.end(), i) - V.begin() + 1;
	    }
    }
    
    // A - class of equiv
    // L - sorted order of suffices
    vector<int> A, L;
    
    int get(int x) { return x >= N ? 0 : A[x]; }
    
    // stable sort elements in L in a by b
    void sort_by(int x) { 
        vector<int> cnt(N + 1); 
        for (int i = 0; i < N; ++i) cnt[get(i + x)]++;
        int sum = 0; 
        for (int i = 0; i <= N; ++i) cnt[i] = (sum += cnt[i], sum - cnt[i]);
        
        vector<int> L2(N);
        for (int i: L) L2[cnt[get(i + x)]++] = i;
        swap(L, L2);
    }
    
    void init(string _str) {
        str = _str; 
        N = str.size();
        
        A.resize(N); 
        for (int i = 0; i < N; ++i) A[i] = str[i]; 
        compress(A);

        L.resize(N); 
        for (int i = 0; i < N; ++i) L[i] = i;
        
        for (int cnt = 1; cnt < N; cnt <<= 1) { 
            sort_by(cnt); 
            sort_by(0);
        
            vector<int> A2(N);
            for (int i = 0; i < N; ++i) {
                if (i == 0) A2[L[i]] = 1;
                else A2[L[i]] = A2[L[i - 1]] +
                    (mp(get(L[i]), get(L[i] + cnt)) != mp(get(L[i - 1]), get(L[i - 1] + cnt)));
            }
            
            swap(A, A2);
        }
    }
    
    vector<int> lcp(vector<int> &inv) { 
        int n = str.size(), h = 0;
        inv.assign(n, 0);
        vector<int> res(n);
        for (int i = 0; i < N; ++i) inv[L[i]] = i;
        for (int i = 0; i < N; ++i) if (inv[i]) {
            int p0 = L[inv[i] - 1];
            while (max(i, p0) + h < N && str[i + h] == str[p0 + h]) h++;
            res[inv[i]] = h;
            if (h) h--;
        }
        return res;
    }
};

void solve() {
    string a = "abacabnaaa";
    cout << "str: " << a << endl;
    SuffixArray sa; sa.init(a);
    vector<int> inv;
    vector<int> lcp = sa.lcp(inv);
    print_vector("lcp", lcp);
    print_vector("inv", inv);
    
    SegTree<int> st(lcp);
    for (int l = 0; l < a.size(); ++l) {
        for (int r = l + 1; r < a.size(); ++r) {
            int a = inv[l];
            int b = inv[r];
            if (a > b) swap(a, b);
            a++; b++;
            cout << l << ' ' << r << ' ' << st.get(a, b - 1) << endl;
        }
    }
}

int main() {
    solve();
    return 0;
}
