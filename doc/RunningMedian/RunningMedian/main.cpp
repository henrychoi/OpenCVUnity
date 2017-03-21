#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

template <typename T> class SortedQ {
    class dValue {
        friend class SortedQ;
    protected:
        dValue *_prev, *_next;
    public:
        T val;
        bool stale;
        inline void init(T v) { stale = false; val = v; }
        inline dValue* next() {
            return (!_next || !_next->stale) ? _next : _next->_next;
        }
        inline dValue* prev() {
            return (!_prev || !_prev->stale) ? _prev : _prev->_prev;
        }
        static int compare(const void* lp, const void* rp) {
            const dValue* *l = (const dValue**)lp;
            const dValue* *r = (const dValue**)rp;
            return (*l)->val - (*r)->val;
        }
        inline void rlink(dValue& right) { _next = &right; right._prev = this; }
        inline void llink(dValue& left) { _prev = &left; left._next = this; }
        inline void replace(dValue* old) {
            _prev = old->_prev; _next = old->_next;
            _prev->_next = _next->_prev = this;
        }
        void rinsert(dValue* newnode);
        void linsert(dValue* newnode);
        inline void dropout() {
            if (_next) _next->_prev = _prev;
            if (_prev) _prev->_next = _next;
        }
    };
public:
    static constexpr unsigned SIZE = 5;
    dValue arr[SIZE+1];
    SortedQ() {
        memset(arr, 0, sizeof(arr)); // mass-set pointers to 0
    }
    ~SortedQ() = default;
    void init();
    inline T& median() { return _median->val; }
    
    //Tail is where we will insert the NEXT sample to, so
    //head should be the slot AFTER the tail (circular).
    inline dValue* head() { return &arr[(_tail + 1) % (SIZE+1)]; }
    T& push(T& f);
    void print();
private:
    unsigned _tail;//In stable state, _tail points to the NEXT slot
    dValue* _median;
};

template <typename T>
void SortedQ<T>::dValue::rinsert(SortedQ<T>::dValue* newnode) {
    dValue* n, *prev;
    for(prev = this, n = next(); n; prev = n, n = n->next()) {
        if (n->val > newnode->val) {
            n->_prev = newnode; //take care of this while I sucked up the branch
            break;
        }
    }
    // Now insert newnode between prev and n
    newnode->_next = n; newnode->_prev = prev;
    prev->_next = newnode; //n->_prev already taken care of above
}
template <typename T>
void SortedQ<T>::dValue::linsert(SortedQ<T>::dValue* newnode) {
    dValue* n, *next;
    for(next = this, n = prev(); n; next = n, n = n->prev()) {
        if (n->val < newnode->val) {
            n->_next = newnode; //take care of this while I sucked up the branch
            break;
        }
    }
    // Now insert newnode between n and next
    newnode->_prev = n; newnode->_next = next;
    next->_prev = newnode; //n->_next already taken care of above
}

template <typename T>
void SortedQ<T>::init() {
    // SIZE number of elements are loaded in q; sort them
    dValue* parr[SIZE];
    for (unsigned i=0; i < SIZE; ++i) {
        parr[i] = &arr[i];
    }
    qsort(parr, SIZE, sizeof(parr[0]), dValue::compare);
    _median = parr[SIZE/2];
    
    for (unsigned i=0; i < SIZE-1; ++i)
        parr[i]->rlink(*parr[i+1]);
    _tail = SIZE;
}

template <typename T>
T& SortedQ<T>::push(T& f) {
    dValue* hole = head(), *tail = &arr[_tail];
    tail->init(f);
    hole->stale = true;//The oldest (head) becomes a stale now
    
    // Find the insertion point for the new value, depending on the decision tree
    if (_median->stale) {
        if (f < _median->_prev->val) { // Left growing; Move _median left
            dValue* prev = _median->_prev;
            prev->linsert(tail);
            _median = prev;
            hole->dropout();
        } else if (_median->_next->val < f) { // Right growing; move median right
            dValue* next = _median->_next;
            next->rinsert(tail);
            _median = next;
            hole->dropout();
        } else { // f replaces median
            tail->replace(_median);
            _median = tail;
        }
    } else { //median is not stale; something to the left or right are stale
        T& med = median();
        if (f >= med) { //f >= _median: stick it to the right
            _median->rinsert(tail);
            if (hole->val < med) { // Hole is to the left; move median right
                _median = _median->_next;
            }
            // What about the equals case?
        } else { //f < _median: stick it to the left
            _median->linsert(tail);
            if (med < hole->val) { // hole is to the right; move median left
                _median = _median->_prev;
            }
            // What about the equals case?
        }
        hole->dropout();
    }
    
    _tail = (_tail + 1) % (SIZE+1); // slot rotates
    return hole->val;
}
template <typename T>
void SortedQ<T>::print() {
    dValue* n, *next;
    for(next = _median, n = _median->_prev; n; next = n, n = n->_prev);
    // Now next points to the first node
    for(n = next; n; n = n->_next) {
        printf((n == _median) ? "%.0f*": "%.0f ", n->val);
    }
}

int main(int argc, const char * argv[]) {
    SortedQ<float> q;
    float f;
    for (unsigned i=0; i < q.SIZE; ++i) {
        f = rand()%4;
        q.arr[i].val = f;
    }
    q.init();
    printf("%.0f [", q.median());
    q.print();
    printf("]\n");
    
    for (unsigned i=0; i < 10; ++i) {
        float oldest = q.push(f = rand()%4);
        printf("%.0f --> %.0f [", f, q.median());
        q.print();
        printf("] --> %.0f\n", oldest);
    }
    
    return 0;
}
