#pragma once
#include <helper/SSharedPtr.hpp>

SNSBEGIN

class SAeroWindow: public SWindow
{
	DEF_SOBJECT(SWindow,L"aero")
public:
	SAeroWindow(void);
	~SAeroWindow(void);

public:
	SOUI_ATTRS_BEGIN()
		ATTR_INT(L"radius",m_radius,TRUE)
	SOUI_ATTRS_END()
protected:
	virtual void OnCommitSurface(IRenderTarget *pRtDest,LPCRECT pRcDest,IRenderTarget *pRtSrc,LPCRECT pRcSrc,BYTE alpha) OVERRIDE;
protected:
	void OnSize(UINT nType, CSize size);
	SOUI_MSG_MAP_BEGIN()
		MSG_WM_SIZE(OnSize)
	SOUI_MSG_MAP_END()

private:
	int m_radius;
	SAutoRefPtr<IRenderTarget> m_midRt;

	class BytePtrDisposer: public PtrDisposer<BYTE>
	{
	public:
		/**
		* Delete pointer ptr.
		* @param ptr pointer to be deleted.
		*/
		virtual void dispose(BYTE *ptr){
			if(ptr)
			{
				delete []ptr;
			}
		}
	};

	
	SSharedPtr<BYTE,BytePtrDisposer> m_buf;
};

SNSEND