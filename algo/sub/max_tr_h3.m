% ハンガリー法による割当問題解（サイズの異なる場合）

function [T] = max_tr_h3(X)
    rows = size(X,2);
    cols = size(X,1);
    T = Hungarian3(X',zeros(rows,cols),cols,rows);
    T = sparse((reshape(T,rows,cols))');
end
