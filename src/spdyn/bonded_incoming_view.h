#pragma once
struct MoveField {
    void *target;
    const void *source;
    int rows,capacity,source_rows,source_capacity,source_offset,bytes;
};
struct MoveView {
    int *counts;
    const int *adds,*exits,*indices;
    int cells,index_capacity,fields;
    MoveField field[12];
};
